
setwd("C:/Users/torn/Documents/Ramfjorden paper")

Data<-read.table("Olavslight_all.txt", header=T)
str(Data)
Dataagg<-aggregate(Data, by=list(Data$Date), FUN=mean)


Data<-read.table("Rampro_weather.txt", header=T)
str(Data)
Data$day[Data$Date == "13.12."]


####################################

#################################################################



abline(h=0, col="blue", lwd=2)
points(c(6,20,62,82,111,139,174,193), c(rep(0,8)), pch=21, bg="red", cex=2)
legend(130,10, c("Temp Ramfjorden", "Temp Tromsø", "normal Temp","sampling time"), lty=c(1,2,1,NA), 
       col=c("black", "red", "grey","black"), pch=c(NA,NA, NA,21), pt.bg=c(NA,NA,NA,"red"))
title(main="(a)",adj=0, cex=0.75)



pdf(file="Rampro_PolNight_Fig2.pdf", height=15, width=15)

par(mar=c(5.1,5.1,4.1,4.1))
par(mfrow=c(2,2))

plot(Data$Temp_RF ~ Data$day, type="l", axes=F, 
     xlab="Date", ylab="Mean daily temperature [°C]", lwd=1.5, ylim=c(-20,11), cex.lab=1.3)
points(Data$Temp_Tr ~ Data$day, type="l", col="red", lwd=0.7, lty=2)
points(Data$normtemp ~ Data$day, type="l", col="grey", lwd=2, lty=1)
axis(2, at=c(-20,-15,-10,-5,0,5,10), las=1, cex.axis=1.3)
axis(1, at=c(11,42,72,103,134,162,193), labels=FALSE)
mtext(c("1.Oct","1.Nov","1.Dec","1.Jan","1.Feb","1.Mar","1.Apr"),1, 
      at=c(11,42,72,103,134,162,193), line=0.5, cex=1, las=2)
abline(h=0, col="blue", lwd=2)
points(c(6,20,62,82,111,139,174,193), c(rep(0,8)), pch=21, bg="red", cex=2)
legend(15,-10, c("Temp Ramfjorden", "Temp Tromsø", "Normal Temp","Sampling time"), lty=c(1,2,1,NA), 
       col=c("black", "red", "grey","black"), pch=c(NA,NA, NA,21), pt.bg=c(NA,NA,NA,"red"), cex=1.3)
title(main="A",adj=0, cex=0.75)

par(mar=c(5.1,5.1,4.1,4.1))
plot(Data$prec ~ Data$day, type="n", axes=F, 
     xlab="Date", ylab="Daily precipitation [mm]", ylim=c(0,25), cex.lab=1.3)
polygon(c(Data$day[-102], rev(Data$day[-102])), c(Data$snow[-102]/10, rep(0,length(Data$snow[-102]))),
        col="lightgrey", border=NA)
polygon(c(Data$day[-102], rev(Data$day[-102])), c(Data$prec[-102], rep(0,length(Data$prec[-102]))),
        col=rgb(0, 0, 0.75, 0.15), border=NA)
points(Data$prec ~ Data$day, type="l", 
      lwd=1.5, col="grey")
axis(2, las=1, cex.axis=1.3)
axis(4, at=c(0,5,10,15,20),labels=FALSE, cex.axis=1.3)
mtext(c("0","50","100","150", "200"),4, line=0.5,cex=1, at=c(0,5,10,15, 20)) 
mtext(c("Snow depth [cm]"), 4, line=1.5,at=c(10))
axis(1, at=c(11,42,72,103,134,162,193), labels=FALSE)
mtext(c("1.Oct","1.Nov","1.Dec","1.Jan","1.Feb","1.Mar","1.Apr"),1, 
      at=c(11,42,72,103,134,162,193), line=0.5, cex=1, las=2)
abline(h=0, col="grey", lwd=2)
points(c(6,20,62,82,111,139,174,193), c(rep(0,8)), pch=21, bg="red", cex=2)
legend(20,24, c("Precipitation Tromsø", "Snow depth Tromsø","Sampling time"), 
       col=c("grey","grey","black"), pch=c(23, 23, 21), pt.bg=c(rgb(0, 0, 0.75, 0.15),"lightgrey","red"), cex=1.3)
title(main="B",adj=0, cex=0.75)

plot(Data$maxwind ~ Data$day, type="n", axes=F, 
     xlab="Date", ylab="Wind speed [m/s]", ylim=c(0,12.8), cex.lab=1.3)
polygon(c(Data$day[-102], rev(Data$day[-102])), c(Data$maxwind[-102], rep(0,length(Data$maxwind[-102]))),
        col=rgb(1, 0, 0, 0.3), border=NA)
polygon(c(Data$day[-c(27,102)], rev(Data$day[-c(27,102)])), c(Data$wind[-c(27,102)], rep(0,length(Data$wind[-c(27,102)]))),
        col="grey", border=NA)
points(Data$wind ~ Data$day, col="black", type="l")
axis(2, las=1, cex.axis=1.3)
axis(1, at=c(11,42,72,103,134,162,193), labels=FALSE)
mtext(c("1.Oct","1.Nov","1.Dec","1.Jan","1.Feb","1.Mar","1.Apr"),1, 
      at=c(11,42,72,103,134,162,193), line=0.5, cex=1, las=2)
abline(h=0, col="grey", lwd=2)
points(c(6,20,62,82,111,139,174,193), c(rep(0,8)), pch=21, bg="red", cex=2)
legend(1,12.5, c("Average wind speed", "Maximum wind speed","Sampling time"), 
       col=c("black","grey","black"), pch=c(23, 23, 21), pt.bg=c("grey",rgb(1, 0, 0, 0.3),"red"), cex=1.3)
title(main="C",adj=0, cex=0.75)



Data2<-Data[ Data$day > 55 & Data$day < 148 & 
               Data$day != "66" & Data$day != "67" & Data$day != "68" &
               Data$day != "78" & Data$day != "79" & Data$day != "80" & Data$day != "81" & Data$day != "82" & Data$day != "83"
             ,]
plot(Data2$PAR ~ Data2$day, type="n", axes=F, 
     xlab="Date", ylab=expression(paste("Average daily Light  PAR  [ ",mu,"E m"^"-2","s"^"-1","]")), cex.lab=1.3)
polygon(c(Data2$day, rev(Data2$day)), c(Data2$PAR, rep(0,length(Data2$PAR))),
        col="gold", border=NA)
points(Data2$PAR ~ Data2$day, col="black", type="l")
axis(2, las=1, cex.axis=1.3)
axis(1, at=c(50,72,103,134,180), labels=FALSE)
mtext(c("1.Dec","1.Jan","1.Feb"),1, 
      at=c(72,103,134), line=0.5, cex=1, las=2)
abline(h=0, col="grey", lwd=2)
points(c(62,82,111,139), c(rep(-0.2,4)), pch=21, bg="red", cex=2)
points(c(62,62,82,82,82,111,111,139), 
       c(2.17,2.54,0.08,0.20,0.13,0.49,0.2,8.5), pch=25, bg="yellow", cex=2)
legend(70,12.5, c("PAR Olavsvern", "PAR at 5m depth","Sampling time"), 
       col=c("black","black","black"), pch=c(23, 25, 21), pt.bg=c("gold","yellow","red"), cex=1.3)
title(main="D",adj=0, cex=0.75)

dev.off()
