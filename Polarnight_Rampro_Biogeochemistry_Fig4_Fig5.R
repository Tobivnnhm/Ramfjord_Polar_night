### load data
getwd()
#C:/Users/torn/Documents"
setwd("C:/Users/torn/Documents/Ramfjorden paper")
Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)
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
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$NOx[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$NOx[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,2,4,6,8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="A",adj=0, cex=0.75)

#### POC/PON at 5m
plot(as.numeric(Data$C.N_ratio_[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("C:N [mol:mol]")), xlim=c(0.5,8.2), ylim=c(0,40),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$C.N_ratio_[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$C.N_ratio_[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,10,20,30,40))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="B",adj=0, cex=0.75)


#### d13C at 5m
plot(as.numeric(Data$d13C[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste(delta^"13","C [","\211","]")), xlim=c(0.5,8.2), ylim=c(-35.5,-25),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$d13C[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$d13C[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(-35,-30,-25))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.75, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="C",adj=0, cex=0.75)

#### NOx at 30m
plot(as.numeric(Data$NOx[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NOx [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$NOx[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$NOx[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,2,4,6,8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="D",adj=0, cex=0.75)


#### POC/PON at 30m
plot(as.numeric(Data$C.N_ratio_[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("C:N [mol:mol]")), xlim=c(0.5,8.2), ylim=c(0,60),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$C.N_ratio_[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$C.N_ratio_[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,20,40,60))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="E",adj=0, cex=0.75)



#### d13C at 30m
plot(as.numeric(Data$d13C[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste(delta^"13","C [","\211","]")), xlim=c(0.5,8.2), ylim=c(-35.5,-25),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$d13C[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$d13C[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(-35,-30,-25))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.75, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
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
     ylab=expression(paste("Chl [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,12.6),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$Chl[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$Chl[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,2,4,6,8,10,12))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="A",adj=0, cex=0.75)


#### POC at 5m
plot(as.numeric(Data$POC[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("POC [",mu,"g L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.8),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$POC[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$POC[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="B",adj=0, cex=0.75)



#### PON at 5m
plot(as.numeric(Data$d15N[Data$station == "st1" & Data$Depth == "5"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "5" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste(delta^"15","N [","\211","]")), xlim=c(0.5,8.2), ylim=c(0,15),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$d15N[Data$station == "st2" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "5" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$d15N[Data$station == "st3" & Data$Depth == "5"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "5" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,5,10,15))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="C",adj=0, cex=0.75)


#### Chl at 30m
plot(as.numeric(Data$Chl[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Chl [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,3.2),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$Chl[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$Chl[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,1,2,3))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="D",adj=0, cex=0.75)



#### POC at 30m
plot(as.numeric(Data$POC[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("POC [",mu,"g L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.8),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$POC[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$POC[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="E",adj=0, cex=0.75)





#### Bacteria at 30m
plot(as.numeric(Data$bacteria[Data$station == "st1"  ])/100000~
       as.numeric(Data$Month2[Data$station == "st1"  ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Bacteria 10"^"5","[cells mL"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$bacteria[Data$station == "st2"  ])/100000~
         as.numeric(Data$Month2[Data$station == "st2"  ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$bacteria[Data$station == "st3" ])/100000~
         as.numeric(Data$Month2[Data$station == "st3"  ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
axis(2, las=1, at=c(0,2,4,6,8))
axis(1, at=c(1,2,3,4,5,6,7,8),labels=FALSE)
mtext(c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7,8), line=0.5, cex=0.75, las=2)
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
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1)
title(main="F",adj=0, cex=0.75)

dev.off()















































#######################################################
### potential supplement


#### Phosphate at 5m
str(Data)
plot(as.numeric(Data$Phosphate[Data$station == "st1" & Data$Depth == "30"])~
     as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Phosphate [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.5),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$Phosphate[Data$station == "st2" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
     pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$Phosphate[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
legend(0.5,0.5, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)

#### NOX at 5m
plot(as.numeric(Data$NOX[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("NOX [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,8),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$NOX[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$NOX[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
legend(0.5,8, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)

#### Chl at 5m
plot(as.numeric(Data$Chl[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("Chl [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,3.5),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$Chl[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$Chl[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
legend(0.5,8, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)

#### POC/PON at 30m
plot(as.numeric(Data$C.N_ratio_[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("C.N_ratio_ [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,60),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$C.N_ratio_[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$C.N_ratio_[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
legend(0.5,8, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)


#### d13C at 30m
plot(as.numeric(Data$d13C[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("d13C [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(-35.5,-25),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$d13C[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$d13C[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
legend(0.5,8, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)


#### POC at 30m
plot(as.numeric(Data$POC[Data$station == "st1" & Data$Depth == "30"])~
       as.numeric(Data$Month2[Data$station == "st1" & Data$Depth == "30" ]-0.1),
     axes=F, xlab="Month", 
     ylab=expression(paste("POC [",mu,"mol L"^"-1","]")), xlim=c(0.5,8.2), ylim=c(0,0.8),
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$POC[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$POC[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
     pch=21, bg="burlywood4",cex=1.02)
points(as.numeric(Data$NH4[Data$station == "st2" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st2" & Data$Depth == "30" ]),
       pch=22, bg="burlywood3", cex=1.02)
points(as.numeric(Data$NH4[Data$station == "st3" & Data$Depth == "30"])~
         as.numeric(Data$Month2[Data$station == "st3" & Data$Depth == "30" ]+0.1),
       pch=24, bg="burlywood1", cex=1.02)
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
legend(0.5,8, c("St1 5m", "St2 5m", "St3 5m"), pch=c(21,22,24), 
       pt.bg =c("burlywood4","burlywood3","burlywood1"), cex=1.2)






