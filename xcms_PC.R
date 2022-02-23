# set working directory
-------------------------------------------------------
rm(list = ls()) ## clears entire workspace, place at top of code
getwd() ## tells you where your are in directory space --> copy and paste

dir <- "path/to/directory/"
dat <- "path/to/data/folder"
res <- "path/to/results/folder/"
script <- "path/to/script/folder"
setwd(file.path(dat))

# load required packages
library(tools)
library(MSnbase)
library(msdata)
library(png)
library(magrittr)
library(xcms)
library(CAMERA)
library(rsm)
library(RColorBrewer)
library("viridis") 
library(ggplot2)
library(pcaMethods)
library("paletteer")
library(reshape2)
library(stringi)
library(scico)
#####################################################################################
######################### Load and prepare data  ####################################
#####################################################################################
setwd(file.path(dir))
meta <- read.csv('Meta_PC_mixotrophy.csv', sep=",", stringsAsFactors = FALSE)
season <- meta$season
meta$season <- factor(meta$season, levels = c("september", "october", "november", "decembre", "january", "february"))
sample_code <- meta$sample_code
mzXMLfiles_folder = "PC/"
setwd(file.path(dat))
mzXMLfiles = list.files(mzXMLfiles_folder, recursive = TRUE, full.names = TRUE)
length(mzXMLfiles)  # number of raw data files

# create xcms xset object from all data file; runs WAY faster with multicore tasking enabled;
files<-mzXMLfiles #input files
pd <- data.frame(sample_name = sub(basename(files), pattern = ".mzXML",
                                   replacement = "", fixed = TRUE),
                 sample_group = season,
                 injection_id = seq(1, 76, by = 1),
                 sample = sample_code,
                 stringsAsFactors = FALSE)
View(pd)

#Define a color palette for your data set
group_color <- scico(6, palette = "batlow")
## group_color <- paletteer_d("LaCroixColoR::PassionFruit", n=6) # or inferno
names(group_color) <- c("september", "october", "november", "decembre", "january", "february")
leg <- c("September", "October", "November", "December", "January", "February")

#Load the raw data as an OnDiskMSnExp using MSnbase package
raw_data <- readMSData(files = files, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk", msLevel. = 1, centroided. = TRUE) 

pData(raw_data)
#get all spectra between 16.55 (895) and 16.75 min (915 sec) (ribitol)
ribitol <- raw_data %>%
  filterRt(rt = c(900, 905)) %>%
  spectra
length(ribitol) # how many spectra
sapply(ribitol, fromFile) # from which files? how many files
plot(ribitol[[989]]) # plot data from last spectrum

raw_data %>%
  filterRt(rt = c(990, 1010)) %>%
  filterMz(mz = c(422.2, 422.4)) %>%
  chromatogram(., aggregationFun = "max") %>%
  plot()

raw_data %>%
  filterFile(. , 36) %>%
  filterRt(rt = c(980, 1050)) %>%
  filterMz(mz = c(422.2, 422.4)) %>% 
  plot(., type = "XIC")
## 2 peaks one at 1000 sec the other at 1030 for #332.2, 217.2, 
## 1 peak at 1000 for 395.3 and 422,4

## Plotting BPC
BPC <- raw_data %>% chromatogram(., aggregationFun = "max")
par(mar=c(5,4,3,9))
plot(BPC,col = group_color[raw_data$sample_group], main = "Base peak Chromatograms \nFA fraction, mixotrophy exp")
legend ("right", inset=c(-0.5,0), xpd = TRUE, legend=leg, fill = group_color)

## Get the total ion current by file
tic <- tic(raw_data, initial = FALSE) #Get tic
from <- fromFile(raw_data) #get sample ID
mydata <-data.frame(tic, from) #combine both vector in a dataframe
pd <- pd %>%
  rownames_to_column(var = "rowname")
mydata <- merge(mydata, pd, by.x = "from", by.y = "rowname")
#Produce plot
ggplot(mydata, aes(x=injection_id, y=tic, group=injection_id, fill=sample_group)) +
  scale_fill_scico_d(6, palette="batlow", name= "Season") +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 2e+06)) +
  labs(title = "TIC per file \nPC fraction, mixotrophy experiment", x="", y = "Intensity") 
#####################################################################################
######## Peak-picking & creation of xcmsSet using IPO parameters ####################
#####################################################################################
# Using XCMS online setting for GC-MS
centW.min_peakwidth = 2 # min peak width in sec
centW.max_peakwidth = 8 # maximum peak width in sec (10 on XCMS online)
centW.ppm = 240 #maximum m/z difference tolerated in ppm, 100 ppm = 0.1 difference
centW.mzdiff = 1 # 0.01 on XCMS online
centW.snthresh = 1 # signal to noise ratio cut off (6 on XCMS online)
centW.prefilter = c(3, 100) # 3 peaks with intensity equal or higher than 100 
centW.noise = 0 #inimum intensity required for centroids to beconsidered in the first analysis step 
centW.fitgauss = TRUE
centW.sleep = 1 #number of seconds to pause between plotting peak finding cycles
centW.mzCenterFun = c("apex") # wMean on XMCS online
centW.verbose.columns = TRUE
centW.integrate = 1 #peak limits are found through descenton the mexican hat filtered data

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

## try with ribitol to see if parameters need to be changed:
# peakwidth: final pick 2 to 8
rib_chr <- chromatogram(raw_data, rt = c(990, 1010), mz = c(422.2, 422.3), 
                        aggregationFun="max")
par(mfrow = c(1, 1), mar = c(4, 4.5, 1, 1))
plot(rib_chr[1,24]) #12, 13, 14, 15, 16, 17, 18, 20, 22, 24
plot(xrib)
chromPeaks(xrib)
# ppm
ribitol <- raw_data %>%
  filterRt(rt = c(995, 1005)) %>%
  filterMz(mz = c(422.2, 422.4))
# Extract intensities and split them by file. This will return a list of lists.
ints_per_file <- split(intensity(ribitol), fromFile(ribitol))
#For each file, sum up the intensities.
ints_sum <- lapply(ints_per_file, function(x) sum(unlist(x)))
# plot sum intentisty to get and idea 
ints_sum %>% 
  unlist(.) %>%
  plot()

#' Extract the ribitol data for one file as a data.frame -> the one with the highest total sum of intenstiy in the selected region
rib_df <- as(filterFile(ribitol, file = which.max(ints_sum)), "data.frame")
#' The difference between m/z values from consecutive scans expressed in ppm
diff(rib_df$mz) * 1e6 / mean(rib_df$mz) 
## set ppm accordingly -> process whole dataset

#Creation of the xset by peak picking
xdata <- findChromPeaks(raw_data, param = cwp) 

## frequency of identified peaks per file along the retention time.
setwd(file.path(res))
png(file="MS_23_PC_chromatographicPeakCounds.png",
    width=600, height=350)
plotChromPeakImage(xdata, main = "Chromatographic peak counts")
dev.off()

eic.ribitol <- chromatogram(xdata, mz = c(422.2, 422.4), rt = c(995, 1005))
chromPeaks(eic.ribitol)
plot(eic.ribitol)

ribitol <- xdata %>%
  filterRt(rt = c(990, 1005)) %>%
  filterMz(mz = c(422.2, 422.4)) %>%
  filterFile(., c(1,6, 12, 76)) %>%
  plot(., type="XIC")

#####################################################################################
######## Peak-picking & creation of xcmsSet using IPO parameters ####################
#####################################################################################

# read in group.density parameter values
density.bw = 2 # xcms online is at 10
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
PC <- "S:/My Libraries/My Library/R/MS23/data/2019_seasons_mzXML/PC/"
#PC <- "/home/dolma/SeaDrive/My Libraries/My Library/R/MS23/data/2019_seasons_mzXML/PC/"
setwd(file.path(PC))
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = density.minfrac, 
                        bw = density.bw,
                        minSamples = density.minsamp,
                        binSize = density.mzwid)
xdata_group <- groupChromPeaks(xdata, param = pdp)

# Step 3: Retention time correction
rtp <- ObiwarpParam()
xdata_RT <- adjustRtime(xdata_group, param = rtp) 
plotAdjustedRtime(xdata_RT)
par(mfrow = c(1,2), mar = c(4, 4.5, 0.9, 0.5))
plot(chromatogram(xdata_RT, mz = c(422.2, 422.4), rt = c(995, 1005),
                  adjustedRtime="FALSE"))
plot(chromatogram(xdata_RT, mz = c(422.2, 422.4), rt = c(995, 1005)))

# Step 4: Group again and change bandwidth as desired
## Prepare dry run first
bpc_ribitol <- chromatogram(xdata_RT, mz = c(422.2, 422.4), rt = c(995, 1005),
                           aggregationFun = "max")
plotChromPeakDensity(bpc_ribitol, param = pdp) # dry run to see if the parameters are good

xdata_group2 <- groupChromPeaks(xdata_RT, param = pdp)
head(featureSummary(xdata_group2))
fmat <- featureValues(xdata_group2, value = "into", method = "maxint")
sum(is.na(fmat)) # missing peaks before filling 

#Step 5: Fill missing peaks
xdata_final <- fillChromPeaks(xdata_group2)
sum(is.na(featureValues(xdata_final))) # missing peaks after rescuing
# Save results in csv file
into_final <- featureValues(xdata_final, value = "into", method = "maxint")
mz_final <- 
rt_final <- 
head(xcms_final)

#Another way
library(SummarizedExperiment)
xcms_final <- quantify(xdata_final, filled = FALSE, method = "sum", value = "into")
# extract only integrated signal from detected peaks (filled = FALSE) and sum up the signals of all
 # chromatographic peaks assigned to a feature in a sample (method = sum).
colData(xcms_final) #column annotation = metadata
rowData(xcms_final) # fearures definition
head(assay(xcms_final, "raw")) # features quantification
# let's add another assay to the object, with quantification after filling in the peaks
assays(xcms_final)$raw_filled <- featureValues(xdata_final, method = "sum",
                                        value = "into", filled = TRUE)
## we now have 2 matrices in the feature object, really cool
assayNames(xcms_final) # to get the names of the matrices.
## possible to add more such as normalized data and such

### Final useful info
#' Overview of the performed processings
processHistory(xdata_final)
#' Access the parameter class for a processing step
processParam(processHistory(xdata_final)[[1]])
## Also save in xcms_final
metadata(xcms_final)

#Save xdata_final to speed up any frurther processing
setwd(file.path(res))
save(xdata_final, file = "MS23_xdata_final_PC_season_laptop.RData")
#load("MS23_xdata_final_PC_season_laptop.RData")
# quick quality control test
#increased the margin at the right side by using the par(mar = …) command.
setwd(file.path(res))
png(file="MS_23_PC_TIC.png",
    width=600, height=350)
par(mar=c(5,4,1,7))
boxplot(featureValues(xdata_final, value="into") +1, 
        col = group_color[xdata_final$sample_group], 
        log="y", las=2)
legend ("right", inset=c(-0.25,0), xpd = TRUE, legend=leg, fill = group_color)
dev.off()

#Step 6: convert XCMSnExp back to xcmsSet (XCMSnExp are not compatible with all function and other packages)
xset <- as(xdata_final, "xcmsSet")
sampnames(xset) <- pData(xdata_final)$sample_name
sampclass(xset) <- pData(xdata_final)$sample_group
xset2 <- fillPeaks(xset) 
save(xset2, file = "MS23_xset2_final_PC_season_laptop.RData")

#Step 7: quality control plots
png(file="MS23_PC_mz_deviation.png", width=600, height = 350)
plotQC(xset2, what="mzdevhist")
dev.off()

png(file="MS23_PC_rt_deviation.png", width=600, height = 350)
plotQC(xset2, what="rtdevhist")
dev.off()

png(file="MS23_PC_mz_deviation_along_mz.png", width=600, height = 350)
plotQC(xset2, what="mzdevmass")
dev.off()

png(file="MS23_PC_mz_deviation_along_rt.png", width=600, height = 350)
plotQC(xset2, what="mzdevtime")
dev.off()

png(file="MS23_PC_mz_deviation_per_sample.png", width=600, height = 350)
plotQC(xset2, what="mzdevsample")
dev.off()

png(file="MS23_PC_rt_deviation_per_sample.png", width=600, height = 350)
plotQC(xset2, what="rtdevsample")
dev.off()

#Step 8 : quality control PCA
sdThresh <- 4.0 ## Filter low-standard deviation rows for plot
data <- log(featureValues(xdata_final))+1

png(filename="MS23_PC_PCA_PC1to3.png", width=600, height=600)
pca.result <- pca(data, nPcs=3)
pcap <- plotPcs(pca.result, type="loadings", col = group_color[xdata_final$sample_group])
legend ("top", inset=c(0.2,0), xpd = TRUE, legend=leg, fill = group_color)
dev.off()

#####################################################################################
###################################### CAMERA #######################################
#####################################################################################
#group peaks into pseudospectra
an <- xsAnnotate(xset2,
                 sample=seq(1,length(sampnames(xset2)))) ## constructor; extracts peak table
## Group peaks in to pseudo-spectra
an <- groupFWHM(an, perfwhm = 3) ## group peaks by retention time
an <- findIsotopes(an, mzabs = 0.01)  # optional but recommended, (# annotate isotopic peaks), #mzabs = allowed m/z error
# findIsotopesWithValidation(object = an, ppm = 10,mzabs = 0.01, intval="intb",maxcharge = 3)
#an <- groupCorr(an,graphMethod="lpc",calcIso = TRUE,calcCiS = TRUE,calcCaS = TRUE, cor_eic_th=0.75) #peak group verification with peakshape correlation
  # For every peak inone pseudospectra a pairwise EIC correlation is performed.  
  #If the correlationvalue between two peaks is higher than the thresholdcoreicthit will stay inthe group, otherwise both are separated.  
  #If the peaks are annotated isotopeions, they will not be divided. 
## Select ions
peaks<-getPeaklist(an, intval="into")
peaks<-peaks[order(peaks$pcgroup),]
peaks[is.na(peaks)]<-0

#####################################################################################
########################## Data normalizazion #######################################
#####################################################################################

#Step 1: Remove samples with a ribitol signal below a given threshold (1E6 or 5E5)
Ribitol<-peaks[peaks$mz > 216 & peaks$mz < 218 & peaks$rt < 1005 & peaks$rt > 995,] #select quantification ion for ribitol 

## Plot Ribitol
ribs<-melt(Ribitol, measure.vars = grep('X', colnames(Ribitol)))
png("MS23_PC_ribitol_intensity_217_5e4.png", width = 600, height = 600)
ggplot(ribs, aes(x=variable, y=value)) + geom_bar(stat='identity', width=0.1) + 
  geom_point() + 
  geom_hline(yintercept=5e4) +
  theme(axis.text.x = element_text( size=8, angle=90))
dev.off()
 
## Remove samples from dataset with ribitol areas less then 5E4
remove<-as.vector(ribs[which(ribs$value < 5E4), 'variable'])
peaks_2<-peaks[,!colnames(peaks) %in% remove]

# Step 2: Normalize data with Ribitol ## stopped here, need to check which sample has double ribitol
subset_rt <- peaks_2[peaks_2$rt < 1005 & peaks_2$rt > 995,] #select right rt windows
# check pcgroup info for ions 217, 307, 319, 332 and 422. Ions of ribitol. in this case pcgroup 5
Ribitol<-peaks_2[peaks_2$mz > 216 & peaks_2$mz < 218 & peaks_2$rt < 1005 & peaks_2$rt > 995,] #select quantification ion for ribitol in the selected samples
r1<-as.vector(t(Ribitol[1,grep('X', colnames(Ribitol))])) # extract ribitol value for each remaining samples
data<-data.frame(peaks_2[grep('X', colnames(peaks_2))]) #keep only column with X in the column name
data.norm<-sweep(data, MARGIN = 2, STATS = r1, FUN = '/')
df_norm<-data.frame(mz=peaks_2$mz, rt=peaks_2$rt, pcgroup=peaks_2$pcgroup,data.norm)
dim(df_norm)

# Step 3: Remove Ions Below 150
df_norm<-df_norm[!df_norm$mz < 150,]
dim(df_norm) # 2933 to 79

# Step 4: Remove Internal standard (Cholestane and/or Ribitol) PC Groups
head(df_norm)
df<-df_norm[!df_norm$pcgroup %in% c(5),] # in our case ribitol is pcgroup 5 and there was not cholestane


#####################################################################################
########################## Ion select for untargeted analysis #######################
#####################################################################################

## Step 1: remove any ions below mz > 150
pksLMz<-df[df$mz > 150,]

## Step 2: Subset by frequency 
pksLMz$count<-rowSums(pksLMz[,grep('X', colnames(pksLMz))] !=0) # count number of times entry appears in dataframe 
pksLMz$count # most mass are present in most of the groups
## Select one peak per pcgroup 
pksFreq<-data.frame(pksLMz %>% group_by(pcgroup) %>% top_n(1,count)) #provides only one entry for most metabolites 
pksFreq

## Step 3: Subset by mean 
pksFreq$means<-rowMeans(pksFreq[,grep('X', colnames(pksFreq))])
pksMean<-data.frame(pksFreq %>% group_by(pcgroup) %>% top_n(1,means))
head(pksMean)
range(pksMean$means)

## Step 4: Remove very low concentrated peaks 
pksMean<-pksMean[pksMean$means > 0.001,]
dim(pksMean)

## Subset peaks dataframe 
pksReduced<-data.frame(mz=pksMean$mz, rt=pksMean$rt, pcgroup=pksMean$pcgroup,pksMean[,grep('X', colnames(pksMean))])
head(pksReduced)

#Save
setwd(file.path(res))
write.csv(pksReduced,file="MS23_result_pksReduced_PC.csv", col.names = TRUE, row.names = FALSE)
