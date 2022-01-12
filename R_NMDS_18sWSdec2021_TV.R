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

getwd()

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
                                      
                                      NMDS1<-metaMDS(t(OTURF)) 
                                      plot(NMDS, disp="sites", type="t")
                                      
                                      ### a 
                                      plot(NMDS1, disp="sites", type="n")
                                      points(NMDS1$points[,1], NMDS1$points[,2],
                                             pch=21,
                                             cex=3,
                                             bg = as.character(ENVRF$color))
                                      legend("topleft", unique(ENVRF$month), pch=21, pt.bg=unique(ENVRF$color), cex=1.4)
                                      ##############################################################
                                      ####### 
                                      
                                      ### Watersamples 5m + 30m
                                      setwd("C:/Users/torn/Documents")
                                  
                                      species <- as.data.frame(read_xlsx("R_PN_Barplot.xlsx"))
                                      row.names(species)<-species[,1]
                                      colnames(species)
                                      species2<-species[,2:27]
                                      species2<-species2[rowMeans(species2) > 0,]
                                
                                      relSpecies <- prop.table(t(species2),1)*100 #make abundancy to persentage
                                      relSpecies #relative abundance in persentage
                                      
                                      ENV2<-data.frame(id=colnames(species2), 
                                                       month=substr(colnames(species2), start=1, stop=3),
                                                       depth=substr(colnames(species2), start=5, stop=7),
                                                       col=as.factor(substr(colnames(species2), start=1, stop=3)),
                                                       pch=as.factor(substr(colnames(species2), start=5, stop=7)))
                                      levels(ENV2$col)<-c("grey45", "darkseagreen4","black","darkseagreen3","gold3","gold","yellow")
                                      levels(ENV2$pch)<-c(22,23,21)
                                      
                                      ######################## WS class plot FARGE #############################
                                      
                                      NMDS <- metaMDS(relSpecies,trymac=50, k = 2)
                                      #NMDS <- metaMDS(relSpecies,trymac=50, k = 2, distance="euclidean")
                                      plot(NMDS, disp="sites", type="t")
                                      
                                      colnames(species)<-c("St3 Sep30", "St1 Oct05","St1 Oct30", "St3 Oct30", "St1 Nov05", "St1 Nov30", "St3 Nov30", "St1 Dec05", "St1 Dec30", "St3 Dec30", "St1 Jan05", "St1 Jan30", 
                                                           "St1 Feb05","St1 Feb30", "St1 Mar30")
                                      colnames(species)
                                      
                                      
                                      plot(NMDS, type="n")
                                      ordihull(NMDS, groups=ENV2$depth, draw="polygon", col=c("brown","gold4","orange"))
                                      points(NMDS$points[,1], NMDS$points[,2],
                                             pch=as.numeric(as.character(ENV2$pch)),
                                             cex=2,
                                             bg = as.character(ENV2$col))
                                      legend("topleft", unique(ENV2$month), pch=21, pt.bg=unique(as.character(ENV2$col)), cex=1.5)
                                      legend("topright", c("Net", "5m", "30m"), pch=unique(as.numeric(as.character(ENV2$pch))), pt.bg=c("brown","gold4","orange"), cex=1.5)

                                          
                                      
                                      
                                      
                                      ###### pdf export
                                      
                                      pdf(file="Polarnight_Rampro_Fig8.pdf", height=8, width=16)
                                      
                                    
                                      par(mfrow=c(1,2))
                                      plot(NMDS1, disp="sites", type="n")
                                      points(NMDS1$points[,1], NMDS1$points[,2],
                                             pch=21,
                                             cex=3,
                                             bg = as.character(ENVRF$color))
                                      legend("topleft", unique(ENVRF$month), pch=21, pt.bg=unique(ENVRF$color), cex=1.4)
                                      title(main="(a)",adj=0, cex=1)
                                      
                                      plot(NMDS, type="n")
                                      ordihull(NMDS, groups=ENV2$depth, draw="polygon", col=c("brown","gold4","orange"))
                                      points(NMDS$points[,1], NMDS$points[,2],
                                             pch=as.numeric(as.character(ENV2$pch)),
                                             cex=2,
                                             bg = as.character(ENV2$col))
                                      legend("topleft", unique(ENV2$month), pch=21, pt.bg=unique(as.character(ENV2$col)), cex=1.5)
                                      legend("topright", c("Net", "5m", "30m"), pch=unique(as.numeric(as.character(ENV2$pch))), pt.bg=c("brown","gold4","orange"), cex=1.5)
                                      title(main="(b)",adj=0, cex=1)
                                      
                                      dev.off()
