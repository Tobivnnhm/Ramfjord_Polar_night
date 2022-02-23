#####################################################################################
########################## Set up environment #######################################
#####################################################################################


rm(list = ls()) ## clears entire workspace, place at top of code
getwd() ## tells you where your are in directory space --> copy and paste

dir <- "path/to/directory/"
dat <- "path/to/data/folder"
res <- "path/to/results/folder/"
script <- "path/to/script/folder"
setwd(file.path(dat))

#Load needed packages
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library("FactoMineR")
library("factoextra")
library("paletteer")
library(data.table)
library(plyr)
library(grid)
#####################################################################################
########################## Load data and metadata ###################################
#####################################################################################

# Load and prepare metadata 
setwd(file.path(dir))
meta <- read.csv('Meta_FA_mixotrophy.csv', sep=",", stringsAsFactors = FALSE)
mzXMLfiles_folder = "FA/"
setwd(file.path(dat))
mzXMLfiles = list.files(mzXMLfiles_folder, recursive = TRUE, full.names = TRUE)
length(mzXMLfiles)  # number of raw data files
nothing <- c(7, 11, 29, 32, 16, 19) # samples with nothing detected
mzXMLfiles[nothing] # print file names
empty <- mzXMLfiles[nothing] %>% gsub("FA/", "", .) %>% gsub(".mzXML", "", .)# remove FA/ from name
meta_final <- meta %>% filter(!sample_name %in% empty)# remove empty files from metadata
sample_code <- meta_final$sample_code
season <- meta_final$season

#Load data after peakpicking
setwd(file.path(res))
df <- read.csv("MS23_result_pksReduced_FA.csv", header=TRUE, stringsAsFactors=FALSE)
head(df)
rownames(df)<-paste('PCGRP',df$pcgroup, 'mz', round(df$mz,2), 'rt', round(df$rt, 1), sep='_') #create new rowname with PCgroup, mz and rt
ion_lookup<-df[,colnames(df) %in% c('mz','rt',"pcgroup")] # ion summary
samples<-colnames(df)[grep('X',colnames(df))] # select samples (all column name with a X)
s<-gsub('X','',samples) # remove the X
str(df)
data <- df[,-c(1:3)] # Subset sample data = remove mz, rt and pcgroup
colnames(data) <- s #replace colnames by new names (without the X)
### Pareto scaling
##### Between the 2 previous methods.
##### The data is mean centered and then divided by the square root of the standard deviation for the variable.
##### Larger variables receive more importance that with autoscaling but less than with mean centering alone.
##### Pareto is very often a good chose for processing mass spectromerty data. Larger peaks will generally have better signal/noise
##### than smaller ones and allowing them to have some extra weight, as compared to autoscaling is helpful.
# Function for pareto scaling
paretoscale <- function(z) {
  colmean <- apply(z, 2, mean) # column means
  colsd <- apply(z, 2, sd)  # column standard deviation
  colsqrtsd <- sqrt(colsd) # sqrt of sd
  rv <- sweep(z, 2, colmean,"-")  # mean center
  rv <- sweep(rv, 2, colsqrtsd, "/")  # dividing by sqrtsd
  return(rv)
}

#Log transform, pareto scale data
data <- data %>%
  t()  %>% #to be feed into pca, r = samples, column = variables
  log1p() %>% #to avoid getting NaN when negative values are present
  paretoscale() %>%
  as.data.frame()
head(data)

###############################################################
# PCA process
###############################################################
# Run PCA
pca <- prcomp(data, center=F, scale=F)
# Create a container called "results" for PCA results
results <- summary(pca)

# Extract PCA results into data frames
scree.data <- as.data.frame(results$importance)  # summary table of explained variance
score.data <- as.data.frame(results$x) #score matrix
loadings.data <- as.data.frame(results$rotation)  # loadings matrix

# Save PCA results to file (we'll use later)
setwd(file.path(res))
write.csv(scree.data, "MS23_FA_pca_scree.csv", row.names = TRUE)
write.csv(score.data, "MS23_FA_pca_scores.csv", row.names = TRUE)
write.csv(loadings.data, "MS23_FA_pca_loadings.csv", row.names = TRUE)

# Make a simple scree plot
percent_1 <- results$importance[2,1:10] * 100
percent_2 <- results$importance[3,1:10] * 100

plot(results$importance[2,1:10], type="b", 
     main="Proportion of Explained Variance",
     xlab="PC", ylab="Proportion of Variance")
text(results$importance[2,1:10], labels = percent_1, pos=4, cex=1) 

# Plot the cumulative proportion of variance
plot(results$importance[3,1:10], type="b", 
     main="Cumulative Proportion of Variance",
     xlab="PC", ylab="Proportion of Variance", ylim=c(0, 1))
text(results$importance[3,1:10], labels = percent_2, pos=4, cex=1) 

# Get variance percentages for first 3 PC's
var1 <- round(scree.data[2,1:3] * 100, 1)
var1
###############################################################
# Score plots
###############################################################
# subset to include only PC1 to PC3 scores
data <- score.data[, c(1:3)]
# look at first few rows
data[1:6,1:3]

#Create dataframe to plot, add metadata
indmeta <- merge(meta_final, data, by.x="sample_name", by.y="row.names")
data <- indmeta
# Make custom theme for ggplot
my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=14, face="bold", vjust=2),
                  panel.background = element_rect(fill = 'gray99',
                                                  colour = "black", 
                                                  size=0.5),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=1.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line(),
                  panel.grid.major = element_line(colour = "gray40", linetype="dotted"),
                  panel.grid.minor = element_line(colour = "gray40", linetype="dashed"),
                  legend.justification=c(1,1),
                  legend.position=c(1,1),
                  legend.title=element_blank(),
                  legend.text = element_text(size = 14))

# check variances for PC1 and PC2
var1
# Set color code
data$season <- factor(data$season, levels = c("september", "october", "november", "decembre", "january", "february"))
group_color <- paletteer_d("LaCroixColoR::PassionFruit", n=6) # or inferno
names(group_color) <- c("september", "october", "november", "decembre", "january", "february")
leg <- c("September", "October", "November", "December", "January", "February")

# draw scores plot
# make plot for PC1 vs. PC2
g1 <- ggplot(data, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "gray50") +
  geom_vline(xintercept = 0, colour = "gray50") +
  geom_point(aes(colour = season), size=4, alpha=0.7) +
  ggtitle("PCA Scores Plot") +
  xlab("PC 1 (38.1%)") +
  ylab("PC 2 (11.4%)") +
  scale_color_manual(values=group_color, labels = leg) +
  my.theme +
  theme(legend.position= "right", legend.title = element_text("Season"))

## label oulier samples to identify them using the pc1 and pc2 results
g1 + geom_label_repel(data=data %>% filter(PC1 < -10 & season == "september"), aes(label=sample_code), show.legend = TRUE, force = 4, label.padding = 0.1) 
head(data)
# save as png file
png(file="MS23_FA_scores.plot_PC12.png", height=2400, width=2800, res=350)
g1
dev.off()

g2 <- ggplot(data, aes(PC1, PC3)) +
  geom_hline(yintercept = 0, colour = "gray50") +
  geom_vline(xintercept = 0, colour = "gray50") +
  geom_point(aes(colour = season), size=4, alpha=0.7) +
  ggtitle("PCA Scores Plot") +
  xlab("PC 1 (38.1%)") +
  ylab("PC 3 (11.4%)") +
  scale_color_manual(values=group_color, labels=leg) +
  my.theme +
  theme(legend.position= "right", legend.title = element_text("Season"))
png(file="MS23_FA_scores.plot_PC13.png", height=2400, width=2800, res=350)
g2
dev.off()

g3 <- ggplot(data, aes(PC2, PC3)) +
  geom_hline(yintercept = 0, colour = "gray50") +
  geom_vline(xintercept = 0, colour = "gray50") +
  geom_point(aes(colour = season), size=4, alpha=0.7) +
  ggtitle("PCA Scores Plot") +
  xlab("PC 2 (11.4%)") +
  ylab("PC 3 (8%)") +
  scale_color_manual(values=group_color, labels=leg) +
  my.theme +
  theme(legend.position= "right", legend.title = element_text("Season"))
png(file="MS23_FA_scores.plot_PC23.png", height=2400, width=2800, res=350)
g3
dev.off()
g1
g2
g3
###############################################################
# Loadings plots
###############################################################
loadings <- loadings.data[,c(1:4)] # only keep PC1 to PC3

#draw loading plots
ggplot(loadings, aes(PC1, PC2)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  labs(title = "Loadings plot", x = "PC1", y = "PC2")

# Find important variables (Loadings)
cutof <- 0.18
cutof_pc1u <- 0.1

plot(loadings$PC1, loadings$PC2)
abline(v=cutof_pc1u, col="red")
abline(v=-cutof, col="red")
abline(h=cutof, col="red")
abline(h=-cutof, col="red")

plot(loadings$PC1, loadings$PC3)
abline(v=cutof_pc1u, col="red")
abline(v=-cutof, col="red")
abline(h=cutof, col="red")
abline(h=-cutof, col="red")

# Use the "ifelse" function to mark high loadings
# PC1 loadings
loadings$pc1.change <-
  ifelse(loadings$PC1 > cutof_pc1u,"UP",
         ifelse(loadings$PC1 < -cutof,"DOWN",
                "none"))
# PC2 loadings
loadings$pc2.change <-
  ifelse(loadings$PC2 > cutof,"UP",
         ifelse(loadings$PC2 < -cutof,"DOWN",
                "none"))
# PC3 loadings
loadings$pc3.change <-
  ifelse(loadings$PC3 > cutof,"UP",
         ifelse(loadings$PC3 < -cutof,"DOWN",
                "none"))


# Number of signficant PC1 loadings
length(which(loadings$pc1.change=="UP"))
length(which(loadings$pc1.change=="DOWN"))
length(which(loadings$pc1.change=="none"))

# Number of signficant PC2 loadings
length(which(loadings$pc2.change=="UP"))
length(which(loadings$pc2.change=="DOWN"))
length(which(loadings$pc2.change=="none"))

# Number of signficant PC3 loadings
length(which(loadings$pc3.change=="UP"))
length(which(loadings$pc3.change=="DOWN"))
length(which(loadings$pc3.change=="none"))

# use the "arrange" function in the plyr package to sort
loadings <- arrange(loadings, pc1.change, pc2.change, pc3.change)
g1 <-
  ggplot(loadings, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "red") +
  geom_vline(xintercept = 0, colour = "red") +
  geom_point(size=2.5, pch=21, color="gray20", bg="khaki1") +
  geom_point(data=subset(loadings, pc1.change=="UP"),
             size=4, pch=21, color="black", bg="#4285F4") +
  geom_point(data=subset(loadings, pc1.change=="DOWN"),
             size=4, pch=22, color="black", bg="#F4B400") +
  geom_point(data=subset(loadings, pc2.change=="UP"),
             size=4, pch=23, color="black", bg="#0F9D58") +
  geom_jitter(data=subset(loadings, pc2.change=="DOWN"),
              size=4, pch=24, color="black", bg="#DB4437") +
  ggtitle("PCA Loadings Plot") +
  my.theme

g2 <-
  ggplot(loadings, aes(PC1, PC3)) +
  geom_hline(yintercept = 0, colour = "gray40") +
  geom_vline(xintercept = 0, colour = "gray40") +
  geom_point(size=2.5, pch=21, color="gray20", bg="khaki1") +
  geom_point(data=subset(loadings, pc1.change=="UP"),
             size=4, pch=21, color="black", bg="#4285F4") +
  geom_point(data=subset(loadings, pc1.change=="DOWN"),
             size=4, pch=22, color="black", bg="#F4B400") +
  geom_point(data=subset(loadings, pc3.change=="UP"),
             size=4, pch=23, color="black", bg="#0F9D58") +
  geom_jitter(data=subset(loadings, pc3.change=="DOWN"),
              size=4, pch=24, color="black", bg="#DB4437") +
  ggtitle("PCA Loadings Plot") +
  my.theme

g3 <-
  ggplot(loadings, aes(PC2, PC3)) +
  geom_hline(yintercept = 0, colour = "gray40") +
  geom_vline(xintercept = 0, colour = "gray40") +
  geom_point(size=2.5, pch=21, color="gray20", bg="khaki1") +
  geom_point(data=subset(loadings, pc2.change=="UP"),
             size=4, pch=21, color="black", bg="#4285F4") +
  geom_point(data=subset(loadings, pc2.change=="DOWN"),
             size=4, pch=22, color="black", bg="#F4B400") +
  geom_point(data=subset(loadings, pc3.change=="UP"),
             size=4, pch=23, color="black", bg="#0F9D58") +
  geom_jitter(data=subset(loadings, pc3.change=="DOWN"),
              size=4, pch=24, color="black", bg="#DB4437") +
  ggtitle("PCA Loadings Plot") +
  my.theme

# add text annotations using the grid package

PC1.pos <- grobTree(textGrob("Positively \n correlated \n with PC1",
                             x=0.90, y=0.15, gp=gpar(col="#4285F4", fontsize=8, fontface="bold")))
PC1.neg <- grobTree(textGrob("Negatively \n correlated \n with PC1",
                             x=0.15, y=0.85, gp=gpar(col="#F4B400", fontsize=8, fontface="bold")))
PC2.pos <- grobTree(textGrob("Positively \n correlated \n with PC2",
                             x=0.90, y=0.85, gp=gpar(col="#0F9D58", fontsize=8, fontface="bold")))
PC2.neg <- grobTree(textGrob("Negatively \n correlated \n with PC2",
                             x=0.17, y=0.15, gp=gpar(col="#DB4437", fontsize=8, fontface="bold")))
PC3a.pos <- grobTree(textGrob("Positively \n correlated \n with PC3",
                              x=0.90, y=0.85, gp=gpar(col="#0F9D58", fontsize=8, fontface="bold")))
PC3a.neg <- grobTree(textGrob("Negatively \n correlated \n with PC3",
                              x=0.17, y=0.15, gp=gpar(col="#DB4437", fontsize=8, fontface="bold")))
PC2b.pos <- grobTree(textGrob("Positively \n correlated \n with PC2",
                              x=0.90, y=0.15, gp=gpar(col="#4285F4", fontsize=8, fontface="bold")))
PC2b.neg <- grobTree(textGrob("Negatively \n correlated \n with PC2",
                              x=0.15, y=0.85, gp=gpar(col="#F4B400", fontsize=8, fontface="bold")))

g1a <-
  g1 +
  annotation_custom(PC1.pos) +
  annotation_custom(PC1.neg) +
  annotation_custom(PC2.pos) +
  annotation_custom(PC2.neg)
png(file="MS23_FA_loadings.plot_PC1PC2.png", height=2400, width=2800, res=350)
g1a
dev.off()
g2a <- 
  g2 +
  annotation_custom(PC1.pos) +
  annotation_custom(PC1.neg) +
  annotation_custom(PC3a.pos) +
  annotation_custom(PC3a.neg)

png(file="MS23_FA_loadings.plot_PC1PC3.png", height=2400, width=2800, res=350)
g2a
dev.off()

g3a <- 
  g3 +
  annotation_custom(PC2b.pos) +
  annotation_custom(PC2b.neg) +
  annotation_custom(PC3a.pos) +
  annotation_custom(PC3a.neg)
png(file="MS23_loadings.plot_PC2PC3.png", height=2400, width=2800, res=350)
g3a
dev.off()

g1a
g2a
g3a
# Subset "significant" loadings (bigger or smaller than cutoff)
loadings.sig <- loadings %>%
  subset( PC1 > cutof_pc1u | PC1 < -cutof |
            PC2 > cutof | PC2 < -cutof |
            PC3 > cutof | PC3 < -cutof) # subset significant loadings

loadings.sig[1:6,1:3]  # look at the first few rows

# sanity check - plot the results 
plot(loadings.sig$PC1, loadings.sig$PC2, main="Loadings \n|0.19| cut off",
                  xlab="PC1", ylab="PC2")
plot(loadings.sig$PC1, loadings.sig$PC3, main="Loadings \n|0.19| cut off",
                  xlab="PC1", ylab="PC3")
plot(loadings.sig$PC2, loadings.sig$PC3, main="Loadings \n|0.19| cut off",
                  xlab="PC2", ylab="PC3")

# Write significant loadings to file for later use.
write.csv(loadings.sig, "MS23_FA_sig_loadings_PC123.csv")

# Merge significant PC1 and PC2 loadings with raw data
# Note: use missing values corrected data
data_sig <- df %>%
  as.data.frame()
sig_loading_names <- loadings.sig$Variable
pca.sig.vars <- data_sig[sig_loading_names,]

# Write the results to file for later use.
write.csv(pca.sig.vars, "MS23_FA_dat_sig_loadings.csv")



