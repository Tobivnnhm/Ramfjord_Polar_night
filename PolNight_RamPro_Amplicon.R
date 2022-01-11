##################################################
### Amplicon sequencing analyses BG-2018 ########
### v1.0 Tobias R Vonnahme ######################
### 17.11.2020 Tobias.vonnahme@uit.no ##########
################################################
### Dependencies: PlotAbund.R function #########
###               taxa.pooler.1.4.R function ###
###               vegan package ###############95*10
################################################



#################################################
### import data #################################
#################################################

OTU <- read.table("OTU_table_16S_TOB1.txt", header=T, row.names=1)
taxonomy <- read.table("Taxonomy_16S_TOB1.txt", header=T, sep="\t", row.names=1)
ENV<-read.table("16S_TOB1_env.txt", header=T)
colnames(OTU)<-ENV$ID2

OTU <- read.table("OTU_18S.txt", header=T, row.names=1)
taxonomy <- read.table("Tax_18S.txt", header=T, sep="\t", row.names=1)
ENV<-read.table("ENV_18S.txt", header=T)
colnames(OTU)<-ENV$id3


###############################################
### subsample 18S for protists################
##############################################

levels(taxonomy$Phylum)
Protisttax<-c("Alveolata_unclass", "Apusozoa_unclass", "Bacillariophyta",
              "Cercozoa","Chlorophyta","Ciliophora","Cryptophyta","Dinoflagellata", "Discosea",
              "Euglenozoa","Foraminifera", "Foraminifera_unclass", "Glaucophyta", "Haptophyta",  "Heliozoa",
              "Ochrophyta","Phaeodarea_unclass","prasinophytes_unclass","Protalveolata","Radiozoa","Rhizaria_unclass",
              "Stramenopiles_unclass","Zygnemophyceae_unclass" )
Tax_protist<-taxonomy[taxonomy$Phylum == Protisttax,]
OTU_protist<-OTU[taxonomy$Phylum == Protisttax,]
taxonomy<-Tax_protist
OTU<-OTU_protist

relOTU <- prop.table(t(OTU),1)*100

colnames(OTU)

#################################################
### subset for different projects ##############
###############################################

### 16S

OTURF<-OTU[taxonomy$Genus != "Chloroplast", ENV$project =="Ram"]
dim(OTURF)
ENVRF<-ENV[ ENV$project =="Ram",]
dim(ENVRF)
colnames(OTURF)<-c("Oct","Nov","Jan","Feb","Mar","Apr")
colnames(OTURF)
#levels(as.factor(colnames(OTUBF2)))<-c("XXX"...)
#colnames(OTUBF2) <- c("XX"...) #option to define the sample names on your own
#levels(as.factor(colnames(OTUBF2)))<-c("XXX"...) # option to define another order of samples on your own

#### 18S
OTURF<-OTU[, ENV$project =="RF_PN"]
dim(OTURF)
ENVRF<-ENV[ ENV$project =="RF_PN",]
dim(ENVRF)
colnames(OTURF)<-c("St1 Oct", "St3 Oct","St 1 Nov", "St1 Dec", "St1 Jan", "St1 Feb","St1 Mar", "St1 Apr")
colnames(OTURF)
#levels(as.factor(colnames(OTUBF2)))<-c("XXX"...)
#colnames(OTUBF2) <- c("XX"...) #option to define the sample names on your own
#levels(as.factor(colnames(OTUBF2)))<-c("XXX"...) # option to define another order of samples on your own

######################################################
## stacked barplots #################################
#####################################################

source("taxa.pooler.1.4.R")
dim(OTURF)
Taxa_pooled <- taxa.pooler(data.frame(OTURF,(taxonomy))) # 18S (8 samples, 7 tax levels, no, no)
dim(OTURF)
Taxa_pooled <- taxa.pooler(data.frame(OTURF,(taxonomy[taxonomy$Genus != "Chloroplast",]))) # 16S 6 samples, 6 levels
Taxa_pooled <- taxa.pooler(data.frame(OTURF[taxonomy$Class == "Planctomycetia",],
                                     (taxonomy[taxonomy$Class == "Planctomycetia" & taxonomy$Genus != "Chloroplast",])))


Phylum_all <- t(Taxa_pooled$Phylum)
relPhylum_all <- prop.table(Phylum_all, 2) * 100
Class_all <- t(Taxa_pooled$Class)
relClass_all <- prop.table(Class_all, 2) * 100
Order_all <- t(Taxa_pooled$Order)
relOrder_all <- prop.table(Order_all, 2) * 100
Family_all <- t(Taxa_pooled$Family)
relFamily_all <- prop.table(Family_all, 2) * 100
Genus_all <- t(Taxa_pooled$Genus)
relGenus_all <- prop.table(Genus_all, 2) * 100
Species_all <- t(Taxa_pooled$Species)
relSpecies_all <- prop.table(Species_all, 2) * 100

source("PlotAbund.R")
PlotAbund(relSpecies_all,6)
PlotAbund(relGenus_all,6)
PlotAbund(relFamily_all,3)
PlotAbund(relOrder_all,3)
PlotAbund(relClass_all,3)
PlotAbund(relPhylum_all,5)

relGenus_all[row.names(relGenus_all) == "Candidatus Pelagibacter"]
relGenus_all[row.names(relGenus_all) == "Euzebyaceae"]

dev.off()

##########################################################
### NMDS plots ###########################################
#########################################################

### 18S
ENVRF$color<-ENVRF$month
ENVRF$color<-c("gold","gold","gold3","grey45","black","darkseagreen4", "darkseagreen3", "darkseagreen2")

## 16S
ENVRF$color<-colnames(OTURF)
ENVRF$color<-c("gold","gold3","black","darkseagreen4", "darkseagreen3", "darkseagreen2")


require(vegan)

NMDS<-metaMDS(t(OTURF), dimensions=1) 
plot(NMDS, disp="sites", type="t")

plot(NMDS, disp="sites", type="n")
points(NMDS$points[,1], NMDS$points[,2],
       pch=21,
       cex=2,
       bg = as.character(ENVRF$color))
legend("topright", colnames(OTURF), pch=21, pt.bg=ENVRF$color)
##############################################################
####### RDA

Data<-read.table("Rampro_median_Data_R.txt", header=T)

ENVRF$ID2
Data16S<-Data[Data$Month != "Sep" & Data$Month != "Dec" & Data$Depth == "30" & Data$station == "st1" & Data$spread == "median",]
ENVRF$NOX <- Data16S$NOx
ENVRF$PO4 <- Data16S$Phosphate
ENVRF$NO2 <- Data16S$Nitrite
ENVRF$NH4 <- Data16S$NH4
ENVRF$Chl <- Data16S$Chl
ENVRF$POC <- Data16S$POC
ENVRF$PN <- Data16S$PN
ENVRF$C.N <- Data16S$C.N_ratio_
ENVRF$d13C <- Data16S$d13C
ENVRF$d15N <- Data16S$d15N



##############################################################
### Options for later ########################################
##############################################################
### ANOSIM: Differ the groups significanly in their community structure?
### Alpha diversity (different indices)
### Beta diversity -(OTU turnover)
### Indicator OTUs/ taxa / Differentially abundant taxa
### predict functions (For 16S)
### other biplots (e.g. PCA)
### relate envrinmental variables (RDA)
### ...
############################################################



Data <-read.table("Rampro_multivar_ENV.txt", header=T)
Data <-read.table("Rampro_multivar_ENV_30.txt", header=T)
Data
dim(Data)

plot(Data$Sal ~ Data$Tpdd)
summary(lm(Data$Sal ~ Data$Tpdd))

Driver<-Data[,16:18]
Driver[,1]<-as.numeric(Driver[,1])
row.names(Driver)<-Data[,1]
str(Driver)


ENV<-(Data[,5:15])
ENV$CN<-as.numeric(ENV$CN)
row.names(ENV)<-Data[,1]
str(ENV)

scale.data.frame <- function(dfr) {
  if (!is.data.frame(dfr)) {stop("dfr must be a data frame")}
  x <- dfr
  cols <- sapply(dfr, is.numeric)
  scaledvars <- scale.default(dfr[, cols]) # otherwise we get a recursive loop
  x[, cols] <- scaledvars
  return(x)
}

zENV<-scale.data.frame(ENV)
str(zENV)


m<-lm(as.matrix(zENV)~as.matrix(Driver)) # sign: PAR X Phaeo, Sal X POC, Sal/PAR X PN |  NOx, PO4, NO3, POC, PN X all, 
summary(m)
anova(m)
plot(m)

m<-aov(as.matrix(zENV)~as.matrix(Driver)) # NOx, PO4, NO3, POC, PN signifcant
summary(m)
plot(m)

m<-manova(as.matrix(zENV)~as.matrix(Driver))
summary(m)

m<-envfit(as.matrix(zENV), as.matrix(Driver))
m
str(m)


m<-rda(as.matrix(zENV)~as.matrix(Driver))
m<-rda(as.matrix(zENV)~Driver[,1] + Driver[,2] + Driver [,3])

summary(m) #RDA1 = 36%, RDA2=9%
plot(m)
anova(m) #Driver are correlated to ENV (p=0.002, F=3.15)
str(anova(m))
m
anova(m,by="term")
anova(m,by="axis") #only axis 1 significant (p=0.001, F=7.38)
(anova.cca(m))


anova_axis<-anova(m,by="axis")
anova_terms<-anova<-anova(m,by="term")
anova_axis
anova_terms

