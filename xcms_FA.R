# set working directory
# --------------- Set up working space ---------------------------------------------------------------------
rm(list = ls()) ## clears entire workspace, place at top of code
getwd() ## tells you where your are in directory space --> copy and paste

dir <- "path/to/directory/"
dat <- "path/to/data/folder"
res <- "path/to/results/folder/"
script <- "path/to/script/folder"
setwd(file.path(dat))

# load required packages
library(tools)
library(xcms)
library(CAMERA)
library(rsm)
library(RColorBrewer)
library("viridis") 
library(ggplot2)
library(pcaMethods)
library("paletteer")
library("scico")
library(tidyverse)
#####################################################################################
######################### Load and prepare data  ####################################
#####################################################################################
#Load metadata
setwd(file.path(dir))
meta <- read.csv('Meta_FA_mixotrophy.csv', sep=",", stringsAsFactors = FALSE)

#Create a season vector to group samples
unique(meta$season)
meta$season <- factor(meta$season, levels = c("september", "october", "november", "decembre", "january", "february"))


# Load data
mzXMLfiles_folder = "FA/"
setwd(file.path(dat))
mzXMLfiles = list.files(mzXMLfiles_folder, recursive = TRUE, full.names = TRUE)
length(mzXMLfiles)  # number of raw data files

# create xcms xset object from all data file; runs WAY faster with multicore tasking enabled;
nothing <- c(7, 11, 29, 32, 16, 19) # samples with nothing detected
mzXMLfiles[nothing] # print file names

files <- mzXMLfiles[-nothing] # remove empty files

empty <- mzXMLfiles[nothing] %>% gsub("FA/", "", .) %>% gsub(".mzXML", "", .)# remove FA/ from name
meta_final <- meta %>% filter(!sample_name %in% empty)# remove empty files from metadata
sample_code <- meta_final$sample_code
season <- meta_final$season
pd <- data.frame(sample_name = sub(basename(files), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 injection_id = seq(1, 56, by = 1),
                 sample = sample_code,
                 sample_group = season,
                 stringsAsFactors = FALSE)
View(pd)

#Define a color palette for your data set
#group_color <- scico(6, palette = "batlow")
group_color <- paletteer_d("LaCroixColoR::PassionFruit", n=6) # or inferno
names(group_color) <- c("september", "october", "november", "decembre", "january", "february")
leg <- c("September", "October", "November", "December", "January", "February")

#Load the raw data as an OnDiskMSnExp using MSnbase package
raw_data <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk", msLevel. = 1, centroided. = TRUE) 
pData(raw_data) #check metadata attached to OnDiskMSnExp file
#####################################################################################
######################### Quick look to raw data ####################################
#####################################################################################

#### Plotting base peak chromatograms ####
BPC <- chromatogram(raw_data, aggregationFun = "max") #Extract base peak chromatogram from the rawdata
# Plotting base peak chromatograms
par(mar=c(5,4,3,9))
plot(BPC,col = group_color[raw_data$sample_group], main = "Base peak Chromatograms \nFA fraction, mixotrophy exp")
legend ("right", inset=c(-0.5,0), xpd = TRUE, legend=leg, fill = group_color)
#setwd(file.path(res))
#ggsave("BPC_FA_mixotrophy.png")
dev.off()
#### Plotting total ion chromatograms ####
TIC <- chromatogram(raw_data, aggregationFun = "sum") #Reads data from the files
# Plotting total ion chromatograms
par(mar=c(5,4,3,9))
plot(TIC,col = group_color[raw_data$sample_group], main = "Total ion chromatograms \nPositive mode")
legend ("right", inset=c(-0.5,0), xpd = TRUE, legend=leg, fill = group_color)
#setwd(file.path(res))
#ggsave("TIC_FA_mixotrophy.png")
dev.off()

#### Create boxplots representing the distribution of TIC (total ion currents) per file ####
#Such plots can be very useful to spot problematic or failing MS runs.
#Get the total ion current by file
tic <- tic(raw_data, initial = FALSE) #Get tic
from <- fromFile(raw_data) #get sample ID
mydata <-data.frame(tic, from) #combine both vector in a dataframe
pd <- pd %>%
  rownames_to_column(var = "rowname")
mydata <- merge(mydata, pd, by.x = "from", by.y = "rowname")

#Produce plot
ggplot(mydata, aes(x=injection_id, y=tic, group=injection_id, fill=sample_group)) +
  scale_fill_paletteer_d("LaCroixColoR::PassionFruit") +
  geom_boxplot() +
  labs(title = "TIC per file \nFA fraction, mixotrophy experiment", x="", y = "Intensity", fill = "Season") +
  coord_cartesian(ylim = c(0, 1e+06)) 
  
#ggsave("TIC_per_file_zoom_FA_mixotrophy.png")

#Produce log2 plot
mydata2 <- mydata %>%
  mutate(log = log2(tic)) %>%
  as.data.frame()
ggplot(mydata2, aes(x=injection_id, y=tic, group=injection_id, fill=sample_group)) + 
  geom_boxplot() +
  labs(title = "log2 intentisty per file \nFA fraction, mixotrophy experiment", x = "", y = expression(log[2]~intensity), fill = "Season") +
  scale_fill_paletteer_d("LaCroixColoR::PassionFruit") +
  coord_cartesian(ylim = c(0, 1e+06))
#ggsave("TIC_per_file_log2.png")

#####################################################################################
######## Peak-picking & creation of xcmsSet using XCMS online parameters ####################
#####################################################################################
# Using XCMS online setting for GC-MS
centW.min_peakwidth = 2
centW.max_peakwidth = 10
centW.ppm = 550 #240
centW.mzdiff = 1
centW.snthresh = 1
centW.prefilter = c(3, 100)
centW.noise = 0
centW.fitgauss = TRUE
centW.sleep = 1
centW.mzCenterFun = c("apex")
centW.verbose.columns = TRUE
centW.integrate = 1

# Step 1: Creation of the xset by peak picking: subset to speed up optimization

cwp <- CentWaveParam(ppm = centW.ppm,
                     peakwidth = c(centW.min_peakwidth, centW.max_peakwidth), 
                     snthresh= centW.snthresh, 
                     prefilter = centW.prefilter,
                     mzCenterFun = centW.mzCenterFun,
                     integrate = centW.integrate,
                     mzdiff = centW.mzdiff, 
                     fitgauss = centW.fitgauss, 
                     noise = centW.noise, 
                     verboseColumns = centW.verbose.columns)

# use internal std to fine tune the parameters in this case lauric acid
raw_data %>%
  filterFile(. , 6) %>% # one sample at the time, travel along injection list to make sure rt does not change too much
  filterRt(rt = c(850, 880)) %>%
  filterMz(mz = c(183.1, 183.4)) %>% 
  plot(., type = "XIC")

par(mfrow = c(1,1))
raw_data %>%
  filterRt(rt = c(850, 880)) %>%
  filterMz(mz = c(183.1, 183.4)) %>%
  chromatogram(., aggregationFun = "max") %>%
  plot()

FA_chr <- chromatogram(raw_data, rt = c(850, 880), mz = c(183.1, 183.4), 
                        aggregationFun="max")

par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
plot(FA_chr[1,24]) #12, 13, 14, 15, 16, 17, 18, 20, 22, 24
xFA <- findChromPeaks(FA_chr, param = cwp)
plot(xFA)
chromPeaks(xFA)

# ppm
IS <- raw_data %>%
  filterRt(rt = c(850, 880)) %>%
  filterMz(mz = c(183.1, 183.4))
# Extract intensities and split them by file. This will return a list of lists.
ints_per_file <- split(intensity(IS), fromFile(IS))
#For each file, sum up the intensities.
ints_sum <- lapply(ints_per_file, function(x) sum(unlist(x)))
# plot sum intentisty to get and idea 
ints_sum %>% 
  unlist(.) %>%
  plot()
# reflect the absence of IS in some samples -> good control
#' Extract the IS data for one file as a data.frame -> the one with the highest total sum of intenstiy in the selected region
IS_df <- as(filterFile(IS, file = which.max(ints_sum)), "data.frame")
IS_df <- as(filterFile(IS, file = 50), "data.frame") # check a few other samples
#' The difference between m/z values from consecutive scans expressed in ppm
diff(IS_df$mz) * 1e6 / mean(IS_df$mz) 
## set ppm accordingly -> process whole dataset

#Creation of the xset by peak picking
xdata <- findChromPeaks(raw_data, param = cwp) 
## frequency of identified peaks per file along the retention time.
plotChromPeakImage(xdata, main = "Chromatographic peak counts")
#ggsave("Chromatographix peak counts FA")

eic.IS <- chromatogram(xdata, mz = c(183.1, 183.4), rt = c(850, 880))
chromPeaks(eic.IS)
plot(eic.IS)

xdata %>%
  filterRt(rt = c(850, 880)) %>%
  filterMz(mz = c(183.1, 183.4)) %>%
  filterFile(., c(1,9, 40, 56)) %>%
  plot(., type="XIC")
#####################################################################################
######## Peak-picking & creation of xcmsSet using IPO parameters ####################
#####################################################################################

# read in group.density parameter values
density.bw = 2
density.minfrac = 0.5
density.minsamp = 1
density.mzwid = 0.25

# read in parameter values for selected retcor method
retcor.meth = "obiwarp"
obiwarp.profStep = 1



# specify some additional settings we wish to keep constant
obiwarp.center = NULL
obiwarp.plottype = "deviation" # "none"
density.sleep = 0
loess.plottype = "mdevden" # none

# Step 2: Peak alignment - initial grouping
FA <- "S:/My Libraries/My Library/R/MS23/data/2019_seasons_mzXML/FA/"
setwd(file.path(FA))
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = density.minfrac, 
                        bw = density.bw,
                        minSamples = density.minsamp,
                        binSize = density.mzwid)
xdata_group <- groupChromPeaks(xdata, param = pdp)

# Step 3: Retention time correction
rtp <- ObiwarpParam()
xdata_RT <- adjustRtime(xdata_group, param = rtp) 

par(mfrow = c(1,1))
plotAdjustedRtime(xdata_RT)

par(mfrow = c(1,2), mar = c(4, 4.5, 0.9, 0.5))
xcms::plot(chromatogram(xdata_RT, mz = c(183.1, 183.4), rt = c(850, 880),
                  adjustedRtime="FALSE"))
xcms::plot(chromatogram(xdata_RT, mz = c(183.1, 183.4), rt = c(850, 880)))

# Step 4: Group again and change bandwidth as desired
xdata_group2 <- groupChromPeaks(xdata_RT, param = pdp)
fmat <- featureValues(xdata_group2, value = "into", method = "maxint")
sum(is.na(fmat)) # missing peaks before filling 
#Step 5: Fill missing peaks
xdata_final <- fillChromPeaks(xdata_group2)
sum(is.na(featureValues(xdata_final))) # missing peaks after rescuing

# quick quality control test
#increased the margin at the right side by using the par(mar = …) command.
par(mfrow = c(1,1), mar=c(5,4,1,8))
boxplot(featureValues(xdata_final, value="into") +1, 
        col = group_color[xdata_group2$sample_group], 
        log="y", las=2)
legend ("right", inset=c(-0.45,0), xpd = TRUE, legend=leg, fill = group_color)

#Save xdata_final
setwd(file.path(res))
save(xdata_final, file = "xdata_final_FA_season.RData")

###############################################################################
################################## CAMERA #####################################
###############################################################################
## Since CAMERA has not yet been ported to XCMSnExp,
## we need to convert to xcmsSet. Note that 
## the conversion only makes sense for somple XCMSnSets, 
## without e.g. MS level filtering (where CAMERA would then 
## extract the wrong )
xset <- as(xdata_final, "xcmsSet")
sampnames(xset) <- pData(xdata_final)$sample_name
sampclass(xset) <- pData(xdata_final)$sample_group
xset2 <- fillPeaks(xset) 
save(xset2, file = "MS23_xset2_final_FA_season.RData")

## Prepare some quality control plot
#setwd(res)
#png(file="MS23_PC_mz_deviation.png", width=600, height = 350)
par(mfrow = c(1,1))
plotQC(xset2, what="mzdevhist")

#dev.off()

#png(file="MS23_PC_rt_deviation.png", width=600, height = 350)
plotQC(xset2, what="rtdevhist")
#dev.off()

#png(file="MS23_PC_mz_deviation_along_mz.png", width=600, height = 350)
par(mar=c(4, 6, 6, 15 ), xpd=TRUE)
plotQC_dmi(xset2, what="mzdevmass")


#png(file="MS23_PC_mz_deviation_along_rt.png", width=600, height = 350)
plotQC_dmi(xset2, what="mzdevtime")
#dev.off()

#png(file="MS23_PC_mz_deviation_per_sample.png", width=600, height = 350)
par(mar=c(8,4,4,4))
plotQC_dmi(xset2, what="mzdevsample")
#dev.off()

#png(file="MS23_PC_rt_deviation_per_sample.png", width=600, height = 350)
plotQC_dmi(xset2, what="rtdevsample")
dev.off()

#Step 8 : quality control PCA
sdThresh <- 4.0 ## Filter low-standard deviation rows for plot
data <- log(featureValues(xdata_final))+1

png(filename="MS23_FA_PCA_PC1to3.png", width=600, height=304)
pca.result <- pca(data, nPcs=3)
plotPcs(pca.result, type="loadings", col = group_color[xdata_final$sample_group])
legend ( "topright", legend=leg, fill = group_color, ncol = 2, xpd = TRUE, bty="n")
dev.off()

#####################################################################################
###################################### CAMERA #######################################
#####################################################################################

## Group peaks in to pseudo-spectra
xset_a <- xsAnnotate(xset2, sample = NA)
xset_f <-groupFWHM(xset_a, perfwhm=3)
an <- findIsotopes(xset_f, mzabs = 0.01)  # optional but recommended, (# annotate isotopic peaks), #mzabs = allowed m/z error
# findIsotopesWithValidation(object = an, ppm = 10,mzabs = 0.01, intval="intb",maxcharge = 3)
#an <- groupCorr(an,graphMethod="lpc",calcIso = TRUE,calcCiS = TRUE,calcCaS = TRUE, cor_eic_th=0.75) #peak group verification with peakshape correlation
# For every peak inone pseudospectra a pairwise EIC correlation is performed.  
#If the correlationvalue between two peaks is higher than the thresholdcoreicthit will stay inthe group, otherwise both are separated.  
#If the peaks are annotated isotopeions, they will not be divided. 
## Select ions
## Select ions
peaks<-getPeaklist(an, intval="into")
peaks<-peaks[order(peaks$pcgroup),]
peaks[is.na(peaks)]<-0

#Step 1: Plot Internal standard
#define quantification ions
rt1 = 850
rt2 = 880
mz1 = 183
mz2 = 184
# select quantification ions
IS <-peaks[peaks$mz > mz1 & peaks$mz < mz2 & peaks$rt < rt2 & peaks$rt > rt1,]  
is <-reshape2::melt(IS, measure.vars = grep('X', colnames(IS)))
#png("MS23_FA_IS_intensity_217_5e4.png", width = 600, height = 600)
ggplot(is, aes(x=variable, y=value)) + geom_bar(stat='identity', width=0.1) + 
  geom_point() + 
  geom_hline(yintercept=1e6) +
  theme(axis.text.x = element_text( size=8, angle=90))
#dev.off()


#Step 2: Remove samples from dataset with ribitol areas less then 5E4
#Usually one should remove sanmples with IS below a given treshold (1e6 or other)
# but in this case the IS was not added to all samples. We will keep all of them.
#remove<-as.vector(is[which(is$value < 1e6), 'variable'])
#peaks_2<-peaks[,!colnames(peaks) %in% remove]
# in our case
peaks_2 <- peaks

#Step 3: indentify pc group of IS
peaks_2[peaks_2$mz > mz1 & peaks_2$mz < mz2,] #check several ions
# 183.2 pcgrp 1
# 214.3 pcgrp 1 
# 199.2 pcgrp 1
# 171.2 pcgrp 1

# Step 3: nomalize data with IS
## also not possible in our case as not all samples have ribitol in them
#is1 <- as.vector(t(IS[1, grep('X', colnames(IS))])) # extract ribitol value for each remaining samples
#data <- data.frame(peaks_2[grep('X', colnames(peaks_2))])  #keep only column with X in the column name
#data.norm<-sweep(data, MARGIN = 2, STATS = is1, FUN = '/')
#df_norm<-data.frame(mz=peaks_2$mz, rt=peaks_2$rt, pcgroup=peaks_2$pcgroup,data.norm)
#dim(df_norm)
df_norm <- peaks_2


# Step 4: Remove Ions Below 150
df_norm<-df_norm[!df_norm$mz < 150,]
dim(df_norm)
head(df_norm)

# Step 5: Remove Internal standard (Cholestane and/or Ribitol) PC Groups
df_norm$rt[df_norm$pcgroup == 1] # have look at the ions on pcgroup 1
df_norm$mz[df_norm$pcgroup == 1]
head(df_norm)
df<-df_norm[!df_norm$pcgroup %in% c(1),] # in our case IS is pcgroup 1

# make sure you removed the right one
df %>% subset(., mz > mz1 & mz < mz2 &
                     rt > rt1 & rt < rt2) %>%
  reshape2::melt(., measure.vars = grep('X', colnames(.))) %>%
  ggplot(., aes(x=variable, y=value)) + geom_bar(stat='identity', width=0.1) + 
  geom_point() + 
  geom_hline(yintercept=1e6) +
  theme(axis.text.x = element_text( size=8, angle=90))

#####################################################################################
########################## Ion select for untargeted analysis #######################
#####################################################################################
pksLMz <- df
## Step 1: in how many samples is a pc group detected, add new column to data set with this info
pksLMz$count<-rowSums(pksLMz[,grep('X', colnames(pksLMz))] !=0) # count number of times entry appears in dataframe 
pksLMz$count # most mass are present in most of the groups

## Step 2: Select one peak per pcgroup -> to simplify plotting, you can always access those data afterwards
pksFreq<-data.frame(pksLMz %>% group_by(pcgroup) %>% top_n(1,count)) #provides only one entry for most metabolites 
pksFreq

## Step 3: Subset by mean 
pksFreq$means<-rowMeans(pksFreq[,grep('X', colnames(pksFreq))])
pksMean<-data.frame(pksFreq %>% group_by(pcgroup) %>% top_n(1,means))
head(pksMean)
range(pksMean$means) # from low abundant ions (below 1 to highily abundant ones 4e5)

## Step 4: Remove very low concentrated peaks 
pksMean<-pksMean[pksMean$means > 0.001,]
dim(pksMean)

## Subset peaks dataframe 
pksReduced<-data.frame(mz=pksMean$mz, rt=pksMean$rt, pcgroup=pksMean$pcgroup,pksMean[,grep('X', colnames(pksMean))])
head(pksReduced)

#Save
setwd(file.path(res))
write.csv(pksReduced,file="MS23_result_pksReduced_FA.csv", col.names = TRUE, row.names = FALSE)
