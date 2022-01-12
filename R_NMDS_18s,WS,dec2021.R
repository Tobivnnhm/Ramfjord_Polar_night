##################################################
### Amplicon sequencing analyses BG-2018 ########
### v1.0 Tobias R Vonnahme ######################
### 17.11.2020 Tobias.vonnahme@uit.no ##########
################################################
### Dependencies: PlotAbund.R function #########
###               taxa.pooler.1.4.R function ###
###               vegan package ###############95*10
################################################
setwd("~/1UiT/artikel/r")

require("vegan")
library(readxl)
require("StressStrength")
require("ggplot2")
require("cba")
library(tidyverse)

#################################################
### import data #################################
#################################################


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
Tax_protist<-taxonomy[taxonomy$Phylum %in% Protisttax,]
OTU_protist<-OTU[taxonomy$Phylum %in% Protisttax,]
taxonomy<-Tax_protist
OTU<-OTU_protist

View(Tax_protist)
unique(Tax_protist$Genus, Order==Syndiniales)

relOTU <- prop.table(t(OTU),1)*100
View(relOTU)
colnames(OTU)
unique(Tax_protist$Class)
Tax_protist[(Tax_protist$Class=="Coscinodiscophyceae"),]
#################################################
### subset for different projects ##############
###############################################

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

                                      ##########################################################
                                      ### NMDS plots ###########################################
                                      #########################################################
                                      
                                      ### 18S
                                      ENVRF$color<-ENVRF$month
                                      ENVRF$color<-c("gold","gold","gold3","grey45","black","darkseagreen4", "darkseagreen3", "darkseagreen2")
                                      
                                      require(vegan)
                                      
                                      NMDS<-metaMDS(t(OTURF)) 
                                      plot(NMDS, disp="sites", type="t")
                                      
                                      plot(NMDS, disp="sites", type="n")
                                      points(NMDS$points[,1], NMDS$points[,2],
                                             pch=21,
                                             cex=2,
                                             bg = as.character(ENVRF$color))
                                      legend("bottomright", colnames(OTURF), pch=21, pt.bg=ENVRF$color)
                                      ##############################################################
                                      ####### 
                                      
                                      ### Watersamples 5m + 30m
                                      
                                  
                                      species <- read.table("D:/UiT/1UiT/artikel/r/PN_NMDS_WS_dec2021.txt", header=T, row.names=1)
                                      
                                      relSpecies <- prop.table(t(species),1)*100 #make abundancy to persentage
                                      relSpecies #relative abundance in persentage
                                      
                                      ######################## WS class plot FARGE #############################
                                      
                                      NMDS <- metaMDS(relSpecies,trymac=50, k = 2)
                                      #NMDS <- metaMDS(relSpecies,trymac=50, k = 2, distance="euclidean")
                                      plot(NMDS, disp="sites", type="t")
                                      
                                      colnames(species)<-c("St3 Sep30", "St1 Oct05","St1 Oct30", "St3 Oct30", "St1 Nov05", "St1 Nov30", "St3 Nov30", "St1 Dec05", "St1 Dec30", "St3 Dec30", "St1 Jan05", "St1 Jan30", 
                                                           "St1 Feb05","St1 Feb30", "St1 Mar30")
                                      colnames(species)
                                      
                                      
                                      plot(NMDS, type="n")
                                      points(NMDS$points[,1], NMDS$points[,2],
                                             pch=21,
                                             cex=2,
                                             bg = as.character(ENVRF$color))
                                      legend("bottomright", colnames(species), pch=21, pt.bg=ENVRF$color)

                                                                            
                                    
                                      
                                      
