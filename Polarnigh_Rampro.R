#############################################################
###Fig. 1####################################################
###CTD interpolation plots###################################
##############################################################

### load data
#setwd("C:/Users/torn/Documents")
cells<-read.csv("protistcells_Rampro.csv", sep=";")

#setwd("D:/Tobias")

Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)
str(Datasum)
sort(unique(Datasum$month))

### subset to Sept - Apr
levels(as.factor(Datasum$month))
Datasum<-Datasum[Datasum$month < 7.5,]

#### Find better interpolation package

#####################################################
### Biogeo data ####################################

Data<-read.table("Rampro_Data_R.txt", header=T)

###############################################################
### Fig 2: NOX, C:N, d13C 5m and 30m #########################
#############################################################

pdf(file="Polarnight_Rampro_Fig4.pdf", height=8, width=12)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,3))

#### NOX at 5m
plot(as.numeric(Data$NOx[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NOx [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$NOx[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$NOx[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,2,4,6,8), cex.axis=1.46)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$NOx[Data$Depth == "5" & Data$Month2 == "1" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "2" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "3" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "4" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "5" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "6" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "7" ], na.rm=T),
          median(Data$NOx[Data$Depth == "5" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(0.5,8, c("St1 5 m", "St2 5 m", "St3 5 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="A",adj=0, cex=0.75)

#### POC/PON at 5m
plot(as.numeric(Data$C.N_ratio_[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("C:N [mol:mol]")), xlim=c(0.5,8.2), ylim=c(0,40),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$C.N_ratio_[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$C.N_ratio_[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,10,20,30,40), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "1" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "2" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "3" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "4" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "5" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "6" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "7" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "5" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(0.5,40, c("St1 5 m", "St2 5 m", "St3 5 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="B",adj=0, cex=0.75)


#### d13C at 5m
plot(as.numeric(Data$d13C[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste(delta^"13","C [","\211","]")), xlim=c(0.5,8.2), ylim=c(-35.5,-25),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$d13C[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$d13C[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(-35,-30,-25), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.75, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$d13C[Data$Depth == "5" & Data$Month2 == "1" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "2" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "3" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "4" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "5" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "6" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "7" ], na.rm=T),
          median(Data$d13C[Data$Depth == "5" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(0.5,-32, c("St1 5 m", "St2 5 m", "St3 5 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="C",adj=0, cex=0.75)

#### NOx at 30m
plot(as.numeric(Data$NOx[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NOx [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$NOx[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$NOx[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,2,4,6,8), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$NOx[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T),
          median(Data$NOx[Data$Depth == "30" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(0.5,8, c("St1 30 m", "St2 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="D",adj=0, cex=0.75)


#### POC/PON at 30m
plot(as.numeric(Data$C.N_ratio_[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("C:N [mol:mol]")), xlim=c(0.5,8.2), ylim=c(0,60),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$C.N_ratio_[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$C.N_ratio_[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,20,40,60), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(0.5,60, c("St1 30 m", "St2 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="E",adj=0, cex=0.75)



#### d13C at 30m
plot(as.numeric(Data$d13C[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste(delta^"13","C [","\211","]")), xlim=c(0.5,8.2), ylim=c(-35.5,-25),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$d13C[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$d13C[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(-35,-30,-25), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.75, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$d13C[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(0.5,-25, c("St1 30 m", "St2 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="F",adj=0, cex=0.75)


dev.off()

















pdf(file="Polarnight_Rampro_Fig5.pdf", height=8, width=12)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,3))

##########################################################################

#### Chl at 5m
plot(as.numeric(Data$Chl[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Chl [",mu,"g L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,12.6),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$Chl[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$Chl[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
points(as.numeric(Data$Phaeo[Data$station == "st1" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]),
       pch="-", col="burlywood4", cex=2.46)
points(as.numeric(Data$Phaeo[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch="-", col="burlywood3", cex=2.46)
points(as.numeric(Data$Phaeo[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch="-", col="burlywood1", cex=2.46)
axis(2, las=1, at=c(0,2,4,6,8,10,12), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$Chl[Data$Depth == "5" & Data$Month2 == "1" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "2" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "3" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "4" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "5" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "6" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "7" ], na.rm=T),
          median(Data$Chl[Data$Depth == "5" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(4,12, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="A",adj=0, cex=0.75)


#### POC at 5m
plot(as.numeric(Data$POC[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("POC [",mu,"g L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.8),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$POC[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$POC[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,0.2,0.4,0.6,0.8), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$POC[Data$Depth == "5" & Data$Month2 == "1" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "2" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "3" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "4" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "5" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "6" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "7" ], na.rm=T),
          median(Data$POC[Data$Depth == "5" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(4,0.8, c("St1 5 m", "St2 5 m", "St3 5 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="B",adj=0, cex=0.75)

colnames(cells)
dim(cells)
cells[10,]
cells5mst1<-cells[10,c(2,3,4,5,6)]
cells30mst1<-cells[10,c(7,8,10,12)]
cells30mst3<-cells[10,c(9,11,13,14,15,16)]
#### d15N-PON at 5m
plot(as.numeric(cells5mst1)~
       as.numeric(c(2,3,4,5,6)-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Protists","[cells mL"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,180),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(cells30mst1)~
         as.numeric(c(1,2,3,4)-0.1),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(cells30mst3)~
         as.numeric(c(2,3,4,5,6,7)-0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,30,60,90,120,150,180), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
legend(5,160, c("St1 5 m", "St1 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="C",adj=0, cex=0.75)


#### Chl at 30m
plot(as.numeric(Data$Chl[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Chl [",mu,"g L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,3.2),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$Chl[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$Chl[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
points(as.numeric(Data$Phaeo[Data$station == "st1" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]),
       pch="-", col="burlywood4", cex=2.46)
points(as.numeric(Data$Phaeo[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch="-", col="burlywood3", cex=2.46)
points(as.numeric(Data$Phaeo[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch="-", col="burlywood1", cex=2.46)
axis(2, las=1, at=c(0,1,2,3), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$Chl[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(4,3, c("St1 30 m", "St2 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="D",adj=0, cex=0.75)



#### POC at 30m
plot(as.numeric(Data$POC[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("POC [",mu,"g L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.8),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$POC[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$POC[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,0.2,0.4,0.6,0.8), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$POC[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(4,0.8, c("St1 30 m", "St2 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="E",adj=0, cex=0.75)





#### Bacteria at 30m
plot(as.numeric(Data$bacteria[Data$station == "st1"  ])/100000~
       as.numeric(Data$Month2[Data$station == "st1"  ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Bacteria 10"^"5","[cells mL"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$bacteria[Data$station == "st2"  ])/100000~
         as.numeric(Data$Month2[Data$station == "st2"  ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$bacteria[Data$station == "st3" ])/100000~
         as.numeric(Data$Month2[Data$station == "st3"  ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,2,4,6,8), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$bacteria[ Data$Month2 == "1" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "2" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "3" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "4" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "5" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "6" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "7" ], na.rm=T)/100000,
          median(Data$bacteria[ Data$Month2 == "8" ], na.rm=T)/100000)
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(4,8, c("St1 30 m", "St2 30 m", "St3 30 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="F",adj=0, cex=0.75)

dev.off()






















##############################################################
#############################################################
####3 Production
pdf(file="Polarnight_Rampro_Fig6.pdf", height=8, width=16)

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,2))

setwd("C:/Users/torn/Documents")
prod<-read.csv("Rampro_production_summary.csv", sep=";")
#### NPP
plot(as.numeric(prod$NPP_mgC_m3_d[c(6,7,8,11,12,13,14)])~
       as.numeric(c(4,4,4,5,5,6,6)-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NPP","[mgC m"^"-3","d"^"-1","]")), xlim=c(0.5,6.5), ylim=c(0,4.5),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(prod$NPP_mgC_m3_d[c(1,2,3,4,5,9,10,15,16)])~
         as.numeric(c(1,1,2,2,4,5,5,6,6)-0.1),
       pch=22, bg="burlywood3", cex=1.46)
abline(h=0)
axis(2, las=1, at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6),labels=FALSE)
mtext(c("Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6), line=0.5, cex=1.4, las=2)
abline(v=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), col="grey")
legend(1.5,3, c("St1 5 m", "st1 30 m"), pch=c(21,22), 
       pt.bg =c("burlywood4","burlywood3"), cex=1.5)
title(main="(a)",adj=0, cex=0.75)

#### Glycine
plot(as.numeric(prod$Ref[1:5])~
       as.numeric(c(1,1,2,2,4)-0.1),
     axes=F, xlab="Month", 
     ylab="At 13C/12C [%]", xlim=c(0.5,4.5), ylim=c(1.07,1.15),
     pch=21, bg="black",cex=1.46, cex.lab=1.4)
abline(h=median(prod$Ref[1:5]))
points(as.numeric(prod$X1C.Gly_atomicexcess[1:5])~
         as.numeric(c(1,1,2,2,4)-0.1),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(prod$X2C.Gly_atomicexcess[1:5])~
         as.numeric(c(1,1,2,2,4)-0.1),
       pch=24, bg="burlywood2", cex=1.46)
points(as.numeric(prod$DIC_atomicexcess[1:5])~
         as.numeric(c(1,1,2,2,4)-0.1),
       pch=21, bg="burlywood4", cex=1.46)
axis(2, las=1, at=c(1.07, 1.08,1.09,1.10,1.11,1.12,1.13,1.14,1.15), cex.axis=1.48)
axis(1, at=c(1,2,3,4),labels=FALSE)
mtext(c("Nov","Dec","Jan","Feb"),1, at=c(1,2,3,4), line=0.5, cex=1.4, las=2)
abline(v=c(0.5,1.5,2.5,3.5,4.5), col="grey")
legend(3.2,1.13, c("Reference","1C-Glycine", "2C-Glycine","DIC"), pch=c(21,22,24,21), 
       pt.bg =c("black","burlywood3","burlywood2","burlywood4"), cex=1.5)
title(main="(b)",adj=0, cex=0.75)

dev.off()


















#######################################################
### potential supplement


#### Phosphate at 5m
str(Data)
plot(as.numeric(Data$Phosphate[Data$station == "st1" & Data$Depth == "30"])~
     as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Phosphate [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.5),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$Phosphate[Data$station == "st2" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
     pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$Phosphate[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,0.1,0.2,0.3,0.4,0.5))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$Phosphate[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,0.5, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)



#### NOX at 5m
plot(as.numeric(Data$NOX[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NOX [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$NOX[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$NOX[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,2,4,6,8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$NOX[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$NOX[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$NOX[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$NOX[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$NOX[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$NOX[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$NOX[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,8, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)




#### Chl at 5m
plot(as.numeric(Data$Chl[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Chl [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,3.5),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$Chl[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$Chl[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,0.5,1,1.5,2,2.5,3,3.5))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$Chl[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$Chl[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,8, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)


#### POC/PON at 30m
plot(as.numeric(Data$C.N_ratio_[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("C.N_ratio_ [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,60),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$C.N_ratio_[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$C.N_ratio_[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,20,40,60))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$C.N_ratio_[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,8, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)


#### d13C at 30m
plot(as.numeric(Data$d13C[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("d13C [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(-35.5,-25),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$d13C[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$d13C[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(-35,-30,-25))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$d13C[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$d13C[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,8, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)





#### POC at 30m
plot(as.numeric(Data$POC[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("POC [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.8),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$POC[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$POC[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$POC[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$POC[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,8, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)






#### NH4 at 5m
plot(as.numeric(Data$NH4[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NH4 [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,6),
     pch=21, bg="burlywood4",cex=1.46)
points(as.numeric(Data$NH4[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$NH4[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,1,2,3,4,5,6))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$NH4[Data$Depth == "30" & Data$Month2 == "1" ], na.rm=T),
          median(Data$NH4[Data$Depth == "30" & Data$Month2 == "2" ], na.rm=T),
          median(Data$NH4[Data$Depth == "30" & Data$Month2 == "3" ], na.rm=T),
          median(Data$NH4[Data$Depth == "30" & Data$Month2 == "4" ], na.rm=T),
          median(Data$NH4[Data$Depth == "30" & Data$Month2 == "5" ], na.rm=T),
          median(Data$NH4[Data$Depth == "30" & Data$Month2 == "6" ], na.rm=T),
          median(Data$NH4[Data$Depth == "30" & Data$Month2 == "7" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
legend(0.5,8, c("St1 5m", "st2 5m", "st3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)

#### d15N-PON at 5m
plot(as.numeric(Data$d15N[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste(delta^"15","N [","\211","]")), xlim=c(0.5,8.2), ylim=c(0,15),
     pch=21, bg="burlywood4",cex=1.46, cex.lab=1.47)
points(as.numeric(Data$d15N[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.46)
points(as.numeric(Data$d15N[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.46)
axis(2, las=1, at=c(0,5,10,15), cex.axis=1.48)
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=1, las=2)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5), col="grey")
medval<-c(median(Data$d15N[Data$Depth == "5" & Data$Month2 == "1" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "2" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "3" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "4" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "5" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "6" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "7" ], na.rm=T),
          median(Data$d15N[Data$Depth == "5" & Data$Month2 == "8" ], na.rm=T))
segments(0.5,medval[1],1.5,medval[1],lwd=2)
segments(1.5,medval[2],2.5,medval[2],lwd=2)
segments(2.5,medval[3],3.5,medval[3],lwd=2)
segments(3.5,medval[4],4.5,medval[4],lwd=2)
segments(4.5,medval[5],5.5,medval[5],lwd=2)
segments(5.5,medval[6],6.5,medval[6],lwd=2)
segments(6.5,medval[7],7.5,medval[7],lwd=2)
segments(7.5,medval[8],8.5,medval[8],lwd=2)
legend(4,15, c("St1 5 m", "St2 5 m", "St3 5 m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.5)
title(main="(c)",adj=0, cex=0.75)


##########################################################################################
############################################################################################
### 
setwd("C:/Users/torn/Documents")
PCenv<-read.csv("Meta_PC_mixotrophy.csv")
FAenv<-read.csv("Meta_FA_mixotrophy.csv")

PC<-read.csv("MS23_result_pksReduced_PC.csv")
FA<-read.csv("MS23_result_pksReduced_FA.csv")

str(PC)
unique(PC$X)
dim(PC)
dim(PCenv)
PCenv

#### PC
dim(PC)
colnames(PC)<-substr(colnames(PC),start=2, stop=13)
PC2<-t(PC[,4:79])
colnames(PC2)<-PC$cgroup
identical(PCenv$sample_name, row.names(PC2))
PCenv$sample_name
PCenv$sample_name2<-row.names(PC2)
is.equal(PCenv$sample_name, row.names(PC2))
PCenv$color<-as.factor(PCenv$season)
levels(PCenv$color)<-c("grey45", "darkseagreen", "black", "gold3","gold","yellow")

require(vegan)

NMDS<-metaMDS(PC2)

NMDSpc<-metaMDS(PC2, dimensions=1) 

plot(NMDSpc, disp="sites", type="n")
points(NMDSpc$points[,1], NMDSpc$points[,2],
       pch=21,
       cex=1.5,
       bg = as.character(PCenv$color))
legend("topright", c("Sept","Oct","Nov","Dec","Jan","Feb"), pt.bg=c("yellow","gold","gold3","grey45","black","darkseagreen"), pch=21, cex=1.5)

#### FA
dim(FA)
if (substr(colnames(FA[4]), start=1, stop=1) == "x"){  colnames(FA)<-substr(colnames(FA),start=2, stop=13) }
dim(PC)
FA2<-t(FA[,4:59])
colnames(FA2)<-FA$cgroup
FA2<-FA2[order(as.numeric(substr(row.names(FA2), start=10, stop=12))),]

FAenv$sample_name
FAenv<-FAenv[order(as.numeric(substr(FAenv$sample_name, start=10, stop=12))),]
FAenv$sample_name2<-row.names(FA2)

identical(FAenv$sample_name, row.names(FA2))
#missing samples in composition sheet (FA2 has less)

FAenv$sample_name
row.names(FA2)
# missing are:20190301_053 -> 31
#             20190301_061 -> 37
#             20190308_92  -> 48
#             20190308_95  -> 51
#             20190308_105 -> 59
#             20190308_110 -> 62
identical(FAenv$sample_name[-c(31,37,48,51,59,62)], row.names(FA2))
FAenv<-FAenv[-c(31,37,48,51,59,62),]
identical(FAenv$sample_name, row.names(FA2)) 

FAenv$color<-as.factor(FAenv$season)
levels(FAenv$color)<-c("grey45", "darkseagreen", "black", "gold3","gold","yellow")

require(vegan)

NMDSfa<-metaMDS(FA2)

NMDSfa<-metaMDS(FA2, dimensions=1) 

plot(NMDS, disp="sites", type="n")
points(NMDSfa$points[,1], NMDSfa$points[,2],
       pch=21,
       cex=1.5,
       bg = as.character(FAenv$color))
legend("topright", c("Sept","Oct","Nov","Dec","Jan","Feb"), pt.bg=c("yellow","gold","gold3","grey45","black","darkseagreen"), pch=21, cex=1.5)

###### combi plot
par(mfrow=c(1,2))
plot(NMDSpc, disp="sites", type="n")
points(NMDSpc$points[,1], NMDSpc$points[,2],
       pch=21,
       cex=1.5,
       bg = as.character(PCenv$color))
legend("topright", c("Sept","Oct","Nov","Dec","Jan","Feb"), pt.bg=c("yellow","gold","gold3","grey45","black","darkseagreen"), pch=21, cex=1.5)
plot(NMDSfa, disp="sites", type="n")
points(NMDSfa$points[,1], NMDSfa$points[,2],
       pch=21,
       cex=1.5,
       bg = as.character(FAenv$color))
legend("topright", c("Sept","Oct","Nov","Dec","Jan","Feb"), pt.bg=c("yellow","gold","gold3","grey45","black","darkseagreen"), pch=21, cex=1.5)

