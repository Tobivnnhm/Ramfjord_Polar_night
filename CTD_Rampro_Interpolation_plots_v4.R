setwd("D:/Tobias")

Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)
#Datasum<-read.table("CTD_BF_July_stations.txt", header=TRUE)
Datasum80<-read.table("Chl_Rampro_1980.txt", header=TRUE)
Datasum08<-read.table("Chl_Rampro_2008.txt", header=TRUE)
Datasum08<-read.table("Temp_Rampro_2008.txt", header=TRUE)

Datasum<-Datasum[Datasum$month < 7.04,]

#Datasum$month <- Datasum$dist

range(Datasum[Datasum$month > 7.04,]$month)

require(castr)

########################
### MLD
colnames(Datasum)
mld1<-rep(0, length(unique(Datasum$month)))
for (i in 1:length(mldts)){
  mld1[i]<-mld(Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 2, 7], 
                Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 2, 1])
}

### Strat Index
stratI1<-rep(0, length(unique(Datasum$month)))
#Datasum<-na.omit(Datasum)
for (i in 1:length(stratI1)){
  stratI1[i]<-stratif(Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 1, 7], 
                     Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 1, 1],
                     min.depth=0:10, max.depths=25:150)
}
stratI2<-rep(0, length(unique(Datasum$month)))
#Datasum<-na.omit(Datasum)
for (i in 1:length(stratI1)){
  stratI2[i]<-stratif(Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 2, 7], 
                      Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 2, 1],
                      min.depth=0:10, max.depths=25:1500)
}
stratI3<-rep(0, length(unique(Datasum$month)))
#Datasum<-na.omit(Datasum)
for (i in 1:length(stratI1)){
  stratI3[i]<-stratif(Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 3, 7], 
                      Datasum[Datasum$month == unique(Datasum$month)[i] & Datasum$station == 3, 1],
                      min.depth=0:10, max.depths=25:150)
}


plot(unique(Datasum$month), stratI1, pch=15, col="black", cex=2)
points(unique(Datasum$month), stratI2, pch=16, col="blue", cex=2)
points(unique(Datasum$month), stratI3, pch=17, col="red", cex=2)

# We have a data.frame containing a time series of 41 CTD casts collected
# at the same station over two years, identified by the date of sampling
head(ctd)
str(ctd)
ctd2<-data.frame(Datasum$month, Datasum$depth, Datasum$temp, Datasum$sal, Datasum$F, Datasum$density, Datasum$station)
colnames(ctd2)<-c("date","depth","temp","sal","chla","sigma", "station")

ctd01<-ctd2[ctd2$station == 1,]
ctd02<-ctd2[ctd2$station == 2,]
ctd03<-ctd2[ctd2$station == 3,]

ctd01$dif<-ctd01$depth
for (i in 2:length(ctd01$depth)){
ctd01$dif[i]<-ctd01$depth[i]-ctd01$depth[i-1]
}
ctd01$dif <- replace(ctd01$dif, which(ctd01$dif < 0.12), NA)
ctd01<-na.omit(ctd01)



# Display the data
library("ggplot2")
ggplot(ctd01) + geom_path(aes(x=temp, y=-depth, group=date), alpha=0.6)
# there are a few spikes in temperature

# Despike each profile
library("dplyr")
ctd_clean <- ctd01 %>% group_by(date) %>%
  mutate(temp=despike(temp, mult=3))
ggplot(ctd_clean) + geom_path(aes(x=temp, y=-depth, group=date), alpha=0.6)
# the spikes are gone! Do the same for all variables
# NB: the despiking parameters should really be tested and adjusted
ctd_clean <- ctd01 %>% group_by(date) %>%
  mutate(
    temp=despike(temp, mult=3),
    sal=despike(sal, mult=3),
    chla=despike(chla, k=2, mult=4),
    # NB: chlarescence is more spiky so be less stringent
    sigma=despike(sigma, mult=3)
  )
ggplot(ctd_clean) + geom_path(aes(x=temp, y=-depth, group=date), alpha=0.6)
ggplot(ctd_clean) + geom_path(aes(x=sal, y=-depth, group=date), alpha=0.6)
ggplot(ctd_clean) + geom_path(aes(x=chla, y=-depth, group=date), alpha=0.6)


# Now compute summary statistics on each despiked profile
stats <- ctd_clean %>% group_by(date) %>%
  summarise(
    thermocline = clined(temp, depth, n.smooth=2, k=2),
    pycnocline = clined(sigma, depth),
    strat_index = stratif(sigma, depth, min.depths=0:5, max.depth=20:30),
    DCM = maxd(chla, depth, n.smooth=2, k=3),
    MLD = mld(sigma, depth, ref.depths=0:5, default.depth=25),
    # it is even possible to use variables computed above to make the
    # following computations adapted to each cast:
    # average tempeature in the mixed layer only
    temp_avg = integrate(temp, depth, from=0, to=MLD, fun=mean),
    # stock of Chl a within 10 m of the DCM
    chla_dcm_stock = integrate(chla, depth, from=DCM-10, to=DCM+10)
  )
# Inspect the results
ggplot(stats) + geom_path(aes(x=date, y=-thermocline))
# thermocline depth is still messy despite the smoothing, probably
# related to sharp changes in this coastal region.










require(MBA)
require(autoimage)




DatasumSub1<-Datasum[Datasum$station==3 & Datasum$month == 2.7,]
DatasumSub2<-Datasum[Datasum$station==3 & Datasum$month == 14.47,]
DatasumSub3<-read.table("CTD_Sept20.txt", header=T)

plot(DatasumSub1$depth ~ DatasumSub1$density, ylim=c(50,0), xlim=c(1020,1027), pch=21, bg="black",cex=1,
     xlab="Density [kg/m3]", ylab="Depth [m]")
points(DatasumSub2$depth ~ DatasumSub2$density, bg="red",pch=21,cex=1)
points(DatasumSub3$depth ~ DatasumSub3$density, bg="cyan",pch=21,cex=1)
abline(h=0, col="grey")

################################################################
### F ##################################################

Datasum<-read.table("All_CTD_Rampro_Chl.txt", header=TRUE)

Datasum$Fcal <- Datasum$F
Datasum$Fcal[Datasum$station == 1 & Datasum$month >1.5]<-Datasum$Fcal[Datasum$station == 1 & Datasum$month >1.5]*4.0585-0.2033
Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month != 7.47]<-Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month != 7.47]*4.0585-0.2033
Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month != 7.47]<-Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month != 7.47]*4.0585-0.2033
Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month == 7.47]<-Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month == 7.47]*0.2782+0.0461
Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month == 7.47]<-Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month == 7.47]*0.2782+0.0461
Datasum$Fcal[Datasum$month < 1.5]<-Datasum$Fcal[Datasum$month < 1.5]

Datasum2<-Datasum[Datasum$station==3 & Datasum$flag != "Castway",]

Dataold<-Datasum80[,c(1,2,3)] # Data2<-Datasum[,c(10,1,4)]
mba.intold <- mba.surf(Dataold, 300, 300, extend=FALSE, n=2, m=2, h=3)$xyz.est # h=4 for sal & light/ h=8 for other par
Dataold2<-Datasum08[Datasum08$Depth < 51,] # Data2<-Datasum[,c(10,1,4)]
mba.intold2 <- mba.surf(Dataold2, 300, 300, extend=FALSE, n=2, m=2, h=3)$xyz.est # h=4 for sal & light/ h=8 for other par


Data2<-Datasum2[,c(10,1,11)] # Data2<-Datasum[,c(10,1,4)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=3)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2 & Datasum$flag != "Castway",]
Data3<-Datasum2[,c(10,1,11)] #Data3<-Datasum2[,c(10,1,4)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=3)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1 & Datasum$flag != "Castway",]
Data3<-Datasum2[,c(10,1,11)] #Data3<-Datasum2[,c(10,1,4)]
mba.int3 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=3)$xyz.est # h=4 for sal & light/ h=8 for other par

maxF<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minF<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

maxF<-max(c(mba.intold2$z,mba.intold$z,mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minF<-min(c(mba.intold2$z,mba.intold$z,mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxF
mba.int$z[2,length(mba.int$z[1,])]<- minF
##Image plot
par(mfrow=c(3,2))


par(mfrow=c(2,2))


#### old data

image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("lightgrey","black","red4","chocolate","chocolate1","gold","darkgreen","green","cyan"))(100), axes=F)
points(Dataold[,1], Dataold[,2], pch=19, col="red", cex=0.1)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)


image(mba.intold2, xaxs="r", yaxs="r", ylim=c(30,0), 
      col=colorRampPalette(c("lightgrey","black","red4","chocolate","chocolate1","gold","darkgreen","green","cyan"))(100), axes=F)
points(Dataold2[,1], Dataold2[,2], pch=19, col="red", cex=0.1)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)
rect(0, 30.5, 10, 12, density = NULL, angle = 45,
     col = "white", border = NULL, lwd =0)




#### NEW data
Datasum2<-Datasum[Datasum$station==3 & Datasum$flag != "Castway",]
Data2<-Datasum2[,c(10,1,11)] # Data2<-Datasum[,c(10,1,4)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

par(mfrow = c(2,2))
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("lightgrey","black","red4","chocolate","chocolate1","gold","darkgreen","green","cyan"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="red", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), las=2, line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)


Datasum2<-Datasum[Datasum$station==2 & Datasum$flag != "Castway",]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,11)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxF
mba.int$z[2,length(mba.int$z[1,])]<- minF
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("lightgrey","black","red4","chocolate","chocolate1","gold","darkgreen","green","cyan"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="red", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), las=2, line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==1 & Datasum$flag != "Castway",]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,11)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxF
mba.int$z[2,length(mba.int$z[1,])]<- minF
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("lightgrey","black","red4","chocolate","chocolate1","gold","darkgreen","green","cyan"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="red", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), las=2, line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

# legend

legend.scale(range(mba.int$z, na.rm=T), col=colorRampPalette(c("lightgrey","black","red4","chocolate","chocolate1","gold","green","green","green","cyan","cyan","cyan"))(100), horizontal = TRUE)
mtext("Chlorophyll a [ug L-1]",1, line=2)








################################################
### F for Chl calib#############################
################################################
Month <- 3.35
station <-1
median(Datasum$F[Datasum$month == Month & Datasum$depth < 5.5 & Datasum$depth > 4.5 & Datasum$station == station])
median(Datasum$F[Datasum$month == Month & Datasum$depth < 30.5 & Datasum$depth > 29.5 & Datasum$station == station])
median(Datasum$F[Datasum$month == Month & Datasum$depth < max(Datasum$depth)+0.5 & Datasum$depth > max(Datasum$depth)+1 & Datasum$station == station])
station <- 2
median(Datasum$F[Datasum$month == Month & Datasum$depth < 5.5 & Datasum$depth > 4.5 & Datasum$station == station])
median(Datasum$F[Datasum$month == Month & Datasum$depth < 30.5 & Datasum$depth > 29.5 & Datasum$station == station])
median(Datasum$F[Datasum$month == Month & Datasum$depth < max(Datasum$depth)+0.5 & Datasum$depth > max(Datasum$depth)+1 & Datasum$station == station])
station <- 3
median(Datasum$F[Datasum$month == Month & Datasum$depth < 5.5 & Datasum$depth > 4.5 & Datasum$station == station])
median(Datasum$F[Datasum$month == Month & Datasum$depth < 30.5 & Datasum$depth > 29.5 & Datasum$station == station])
median(Datasum$F[Datasum$month == Month & Datasum$depth < max(Datasum$depth)+0.5 & Datasum$depth > max(Datasum$depth)+1 & Datasum$station == station])

Datasum$Fcal <- Datasum$F
Datasum$Fcal[Datasum$station == 1]<-Datasum$Fcal[Datasum$station == 1]*4.2772-0.2049
Datasum$Fcal[Datasum$station == 2 & Datasum$month != 7.47]<-Datasum$Fcal[Datasum$station == 2 & Datasum$month != 7.47]*4.2772-0.2049
Datasum$Fcal[Datasum$station == 3 & Datasum$month != 7.47]<-Datasum$Fcal[Datasum$station == 3 & Datasum$month != 7.47]*4.2772-0.2049
Datasum$Fcal[Datasum$station == 2 & Datasum$month == 7.47]<-Datasum$Fcal[Datasum$station == 2 & Datasum$month == 7.47]*0.2782+0.0461
Datasum$Fcal[Datasum$station == 3 & Datasum$month == 7.47]<-Datasum$Fcal[Datasum$station == 3 & Datasum$month == 7.47]*4.2782+0.0461


################################################################
### PAR ##################################################
Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)
Datasum2<-Datasum[Datasum$station==3 & Datasum$flag != "Castway",]
Data2<-Datasum2[,c(10,1,5)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2 & Datasum$flag != "Castway",]
Data3<-Datasum2[,c(10,1,5)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1 & Datasum$flag != "Castway",]
Data3<-Datasum2[,c(10,1,5)]
mba.int3 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

maxpar<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minpar<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxpar
mba.int$z[2,length(mba.int$z[1,])]<- minpar

##Image plot
par(mfrow=c(2,2))

image(mba.int, xaxs="r", yaxs="r", ylim=c(120,0), 
      col=colorRampPalette(c("black","red","goldenrod1","orange","yellow",rep("lightyellow",2),rep("white",2)))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Dec","Jan","Feb","Mar","Apr"),1, at=c(3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==2 & Datasum$flag != "Castway",]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,5)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxpar
mba.int$z[2,length(mba.int$z[1,])]<- minpar
image(mba.int, xaxs="r", yaxs="r", ylim=c(80,0), 
      col=colorRampPalette(c("black","red","goldenrod1","orange","yellow",rep("lightyellow",2),rep("white",2)))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Dec","Jan","Feb","Mar","Apr"),1, at=c(3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==3 & Datasum$flag != "Castway" ,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,5)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxpar
mba.int$z[2,length(mba.int$z[1,])]<- minpar
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("black","red","goldenrod1","orange","yellow",rep("lightyellow",2),rep("white",2)))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Dec","Jan","Feb","Mar","Apr"),1, at=c(3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

# legend
require(autoimage)

legend.scale(range(mba.int$z, na.rm=T), 
             col=colorRampPalette(c("black","red","goldenrod1","orange","yellow",rep("lightyellow",2),rep("white",2)))(100),
             horizontal = TRUE)
mtext("PAR [uE m-2 s-1]",1, line=2)





aggregate(Datasum[Datasum$station == 1 & Datasum$depth > 29.5 & Datasum$depth <30.1, c(5,10)],
          by=list(Datasum[Datasum$station == 1 & Datasum$depth > 29.5 & Datasum$depth <30.1, c(5,10)]$month), 
          FUN=mean)
aggregate(Datasum[Datasum$station == 2 & Datasum$depth > 29.5 & Datasum$depth <30.1, c(5,10)],
          by=list(Datasum[Datasum$station == 2 & Datasum$depth > 29.5 & Datasum$depth <30.1, c(5,10)]$month), 
          FUN=mean)
aggregate(Datasum[Datasum$station == 3 & Datasum$depth > 29.5 & Datasum$depth <30.1, c(5,10)],
          by=list(Datasum[Datasum$station == 3 & Datasum$depth > 29.5 & Datasum$depth <30.1, c(5,10)]$month), 
          FUN=mean)


aggregate(Datasum[Datasum$station == 1 & Datasum$depth > 4.5 & Datasum$depth <5.1, c(5,10)],
          by=list(Datasum[Datasum$station == 1 & Datasum$depth > 4.5 & Datasum$depth <5.1, c(5,10)]$month), 
          FUN=mean)


################################################################
### Density ##################################################
Datasum2<-Datasum[Datasum$station==3 & Datasum$month != 7.47,]
Data2<-Datasum2[,c(10,1,7)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2 & Datasum$month != 7.47,]
Data3<-Datasum2[,c(10,1,7)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1 & Datasum$month != 7.47,]
Data3<-Datasum2[,c(10,1,7)]
mba.int3 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

maxS<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minS<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS

##Image plot
par(mfrow=c(2,2))

image(mba.int, xaxs="r", yaxs="r", ylim=c(120,0), 
      col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==2 & Datasum$month != 7.47,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,7)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", ylim=c(80,0), 
      col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==1 & Datasum$month != 7.47,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,7)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

# legend

legend.scale(range(mba.int$z, na.rm=T), 
             col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100),
             horizontal = TRUE)
mtext("Density",1, line=2)

plot(Datasum$depth[Datasum$station == "1" & Datasum$month == "7.03"]~
       Datasum$density[Datasum$station == "1" &  Datasum$month == "7.03"],
     type="l", ylim=c(40,0))


################################################################
### Salinity ##################################################
Datasum2<-Datasum[Datasum$station==3 & Datasum$flag != "Castway",]
Data2<-Datasum2[,c(10,1,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2,]
Data3<-Datasum2[,c(10,1,2)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1,]
Data3<-Datasum2[,c(10,1,2)]
mba.int3 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

maxS<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minS<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS

##Image plot
par(mfrow=c(2,2))

image(mba.int, xaxs="r", yaxs="r", ylim=c(120,0), 
      col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==2,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", ylim=c(80,0), 
      col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==1,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

# legend

legend.scale(range(mba.int$z, na.rm=T), 
             col=colorRampPalette(c(rep("cyan",3),rep("blue",2),rep("yellow",3),"red","brown","darkred","black"))(100),
             horizontal = TRUE)
mtext("Salinity",1, line=2)


################################################################
### Salinity subset##################################################
Datasum2<-Datasum[Datasum$station==3 & Datasum$sal >31,]
Data2<-Datasum2[,c(10,1,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2 & Datasum$sal >31,]
Data3<-Datasum2[,c(10,1,2)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1 & Datasum$sal >31,]
Data3<-Datasum2[,c(10,1,2)]
mba.int3 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

maxS<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minS<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS


Datasum3<-Datasum[ Datasum$sal >31.5,]
Datasum2<-Datasum[Datasum$station==3 & Datasum$sal >31,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,2)]


mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
##Image plot
par(mfrow=c(2,2))

image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==2 & Datasum$sal >31,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", ylim=c(80,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station== 1 & Datasum$sal >31,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

# legend
legend.scale(range(mba.int$z, na.rm=T), 
             col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100),
             horizontal = TRUE)
mtext("Salinity [PSU]",1, line=2)



#################################################################
###### Temperature #############################################


Dataold2<-Datasum08 # Data2<-Datasum[,c(10,1,4)]
mba.intold2 <- mba.surf(Dataold2, 300, 300, extend=FALSE, n=2, m=2, h=3)$xyz.est # h=4 for sal & light/ h=8 for other par








Datasum2<-Datasum[Datasum$station==3,]
Data2<-Datasum2[,c(10,1,3)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2,]
Data3<-Datasum2[,c(10,1,3)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1,]
Data3<-Datasum2[,c(10,1,3)]
mba.int3 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

maxT<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minT<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

maxT<-max(c(mba.intold2$z,mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minT<-min(c(mba.intold2$z,mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxT
mba.int$z[2,length(mba.int$z[1,])]<- minT


image(mba.intold2, xaxs="r", yaxs="r", ylim=c(15,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","orange","red"))(90), axes=F)
points(Dataold2[,1], Dataold2[,2], pch=19, col="black", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)






##Image plot
par(mfrow=c(2,2))
Datasum2<-Datasum[Datasum$station==3,]
Data2<-Datasum2[,c(10,1,3)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
image(mba.int, xaxs="r", yaxs="r", ylim=c(120,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","orange","red"))(90), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="black", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==2,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,3)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxT
mba.int$z[2,length(mba.int$z[1,])]<- minT
image(mba.int, xaxs="r", yaxs="r", ylim=c(80,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","orange","red"))(90), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="black", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==1,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,1,3)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxT
mba.int$z[2,length(mba.int$z[1,])]<- minT
image(mba.int, xaxs="r", yaxs="r", ylim=c(50,0), 
      col=colorRampPalette(c("white","cyan","blue","yellow","orange","red"))(90), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="black", cex=0.03)
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sept","Oct"),1, at=c(1,2,3,4,5,6,7,8,9,10,11,12,13), line=0.5, las=2)
mtext("Depth [m]", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

# legend

legend.scale(range(mba.int$z, na.rm=T), col = colorRampPalette(c("white","cyan","blue","yellow","orange","red"))(90), horizontal = TRUE)
mtext("Temp [C]",1, line=2)



#########################################################
### TS plots

Datasum2<-Datasum[Datasum$station==3 & Datasum$sal >31.5,]
Data2<-Datasum2[,c(10,3,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==2 & Datasum$sal >31.5,]
Data3<-Datasum2[,c(10,3,2)]
mba.int2 <- mba.surf(Data3, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par
Datasum2<-Datasum[Datasum$station==1 & Datasum$sal >31.5,]
Data4<-Datasum2[,c(10,3,2)]
mba.int3 <- mba.surf(Data4, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est # h=4 for sal & light/ h=8 for other par

maxS<-max(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)
minS<-min(c(mba.int$z,mba.int2$z,mba.int3$z),na.rm=T)

rm("mba.int2")
rm("mba.int3")
rm("Datasum2")
rm("Data3")

mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS


##Image plot
par(mfrow=c(2,2))

image(mba.int, xaxs="r", yaxs="r", 
      col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19,
       col="grey", cex=0.03)
polygon(c(min(Data2[,2]),aggregate(Data2[,2], by=list(Data2[,1]), FUN=min)[,1], max(Data2[,2])),
           c(0,aggregate(Data2[,2], by=list(Data2[,1]), FUN=min)[,2],0),
          col="White",border="white")
polygon(c(min(Data2[,2]),aggregate(Data2[,2], by=list(Data2[,1]), FUN=max)[,1], max(Data2[,2])),
        c(12,aggregate(Data2[,2], by=list(Data2[,1]), FUN=max)[,2],12),
        col="White", border="white")
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Temperature", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station==2 & Datasum$sal >31.5,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,3,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r",
      col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
polygon(c(min(Data2[,2]),aggregate(Data2[,2], by=list(Data2[,1]), FUN=min)[,1], max(Data2[,2])),
        c(0,aggregate(Data2[,2], by=list(Data2[,1]), FUN=min)[,2],0),
        col="White",border="white")
polygon(c(min(Data2[,2]),aggregate(Data2[,2], by=list(Data2[,1]), FUN=max)[,1], max(Data2[,2])),
        c(12,aggregate(Data2[,2], by=list(Data2[,1]), FUN=max)[,2],12),
        col="White", border="white")
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Temperature", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

Datasum2<-Datasum[Datasum$station== 1 & Datasum$sal >31.5,]
colnames(Datasum2)
Data2<-Datasum2[,c(10,3,2)]
mba.int <- mba.surf(Data2, 300, 300, extend=FALSE, n=2, m=2, h=4)$xyz.est
mba.int$z[1,length(mba.int$z[1,])]<- maxS
mba.int$z[2,length(mba.int$z[1,])]<- minS
image(mba.int, xaxs="r", yaxs="r", 
      col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100), axes=F)
points(Data2[,1], Data2[,2], pch=19, col="grey", cex=0.03)
polygon(c(min(Data2[,2]),aggregate(Data2[,2], by=list(Data2[,1]), FUN=min)[,1], max(Data2[,2])),
        c(0,aggregate(Data2[,2], by=list(Data2[,1]), FUN=min)[,2],0),
        col="White",border="white")
polygon(c(min(Data2[,2]),aggregate(Data2[,2], by=list(Data2[,1]), FUN=max)[,1], max(Data2[,2])),
        c(12,aggregate(Data2[,2], by=list(Data2[,1]), FUN=max)[,2],12),
        col="White", border="white")
axis(2)
axis(1, labels=FALSE)
mtext(c("Oct","Nov","Dec","Jan","Feb","Mar","Apr"),1, at=c(1,2,3,4,5,6,7), line=0.5)
mtext("Temperature", 2, line=2)
mtext("Month", 1, line=2)
abline(h=0)

legend.scale(range(mba.int$z, na.rm=T), 
             col=colorRampPalette(c("white","cyan","blue","yellow","red","brown","darkred","black"))(100),
             horizontal = TRUE)
mtext("Salinity [PSU]",1, line=2)


#### v2
par(mfrow=c(2,3))
Datasum2<-Datasum[Datasum$station==3 & Datasum$sal >31.5,]
colnames(Datasum2)
Datasum2$col= as.factor(Datasum2$month)
levels(Datasum2$col)<- colorRampPalette(c("red","black","lightgrey","cyan"))(length(levels(Datasum2$col)))
plot(Datasum2$temp ~ Datasum2$sal, col=as.character(Datasum2$col), pch=19, xlab="Salinity [PSU]", ylab="Temperature [C]", main="St 3")

Datasum2<-Datasum[Datasum$station==2 & Datasum$sal >31.5,]
colnames(Datasum2)
Datasum2$col= as.factor(Datasum2$month)
levels(Datasum2$col)<- colorRampPalette(c("red","black","lightgrey","cyan"))(length(levels(Datasum2$col)))
plot(Datasum2$temp ~ Datasum2$sal, col=as.character(Datasum2$col), pch=19, xlab="Salinity [PSU]", ylab="Temperature [C]", main="St 2")


Datasum2<-Datasum[Datasum$station==1 & Datasum$sal >31.5,]
colnames(Datasum2)
Datasum2$col= as.factor(Datasum2$month)
levels(Datasum2$col)<- colorRampPalette(c("red","black","lightgrey","cyan"))(length(levels(Datasum2$col)))
plot(Datasum2$temp ~ Datasum2$sal, col=as.character(Datasum2$col), pch=19, xlab="Salinity [PSU]", ylab="Temperature [C]", main="St 1")

##

Datasum2<-Datasum[Datasum$station==3 & Datasum$sal >31.5,]
colnames(Datasum2)
Datasum2$col= as.factor(Datasum2$depth)
levels(Datasum2$col)<- colorRampPalette(c("white","black"))(length(levels(Datasum2$col)))
plot(Datasum2$temp ~ Datasum2$sal, col=as.character(Datasum2$col), pch=19, xlab="Salinity [PSU]", ylab="Temperature [C]", main="St 3")

Datasum2<-Datasum[Datasum$station==2 & Datasum$sal >31.5,]
colnames(Datasum2)
Datasum2$col= as.factor(Datasum2$depth)
levels(Datasum2$col)<- colorRampPalette(c("white","black"))(length(levels(Datasum2$col)))
plot(Datasum2$temp ~ Datasum2$sal, col=as.character(Datasum2$col), pch=19, xlab="Salinity [PSU]", ylab="Temperature [C]", main="St 2")


Datasum2<-Datasum[Datasum$station==1 & Datasum$sal >31.5,]
colnames(Datasum2)
Datasum2$col= as.factor(Datasum2$depth)
levels(Datasum2$col)<- colorRampPalette(c("white","black"))(length(levels(Datasum2$col)))
plot(Datasum2$temp ~ Datasum2$sal, col=as.character(Datasum2$col), pch=19, xlab="Salinity [PSU]", ylab="Temperature [C]", main="St 2")



#############################
Data<-read.csv("CC1747003_20190314_100333.csv", header=F, skip=1)
colnames(Data)<-c("depth","temp","sal")


plot(Data$sal, Data$depth, type = "l", ylim=c(max(Data$depth),0),  
     
     xlim=c(31,max(Data$sal)),
     
     axes= FALSE, xlab="",ylab="depth [m]")

axis(2,las=1)

axis(1,31:max(Data$sal),line=1,col="black",col.ticks="black",col.axis="black")

mtext("Salinity",1,line=1,at=31,col="black")

par(new=TRUE)

plot(Data$temp, Data$depth, type = "l", ylim=c(max(Data$depth),0), xlim=c(0,max(Data$tem)),
     
     col="blue", axes=FALSE, xlab="",ylab="")

axis(1,0:max(Data$temp), line=3, col="blue", col.ticks="blue",col.axis="blue")

mtext("Temperature",1,line=3,at=0,col="blue")

##############################################################
#Iriarte, A., Villate, F., Uriarte, I., Alberdi, L., & Intxausti, L. (2015). Dissolved oxygen in a temperate estuary: the influence of hydro-climatic factors and eutrophication at seasonal and inter-annual time scales. Estuaries and Coasts, 38(3), 1000-1015.

#1

Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)
levels(as.factor(Datasum$month))

for (i in levels(as.factor(Datasum$month))) {
Subdata<-Datasum[Datasum$station==1 & Datasum$month == i,]
Subdata2<-Datasum[Datasum$station==2 & Datasum$month == i,]
Subdata3<-Datasum[Datasum$station==3 & Datasum$month == i,]

st1<-(Subdata[Subdata$depth == max(Subdata$depth),]$density - 
            Subdata[Subdata$depth == min(Subdata$depth),]$density)/
              mean(Subdata$density)
st2<-(Subdata2[Subdata2$depth == max(Subdata2$depth),]$density - 
    Subdata2[Subdata2$depth == min(Subdata2$depth),]$density)/
  mean(Subdata2$density)
st3<-(Subdata3[Subdata3$depth == max(Subdata3$depth),]$density - 
    Subdata3[Subdata3$depth == min(Subdata3$depth),]$density)/
  mean(Subdata3$density)
st1
st2
st3}

Subdata$dens

######################
### Stratification index 1
##########################


k=19
i=as.numeric(levels(as.factor(Datasum$month)))
Subdata<-Datasum[Datasum$station==1 & Datasum$month == i[k],]
Subdata2<-Datasum[Datasum$station==2 & Datasum$month == i[k],]
Subdata3<-Datasum[Datasum$station==3 & Datasum$month == i[k],]

st1<-(Subdata[Subdata$depth == max(Subdata$depth),]$density - 
        Subdata[Subdata$depth == min(Subdata$depth),]$density)/
  mean(Subdata$density)
st2<-(Subdata2[Subdata2$depth == max(Subdata2$depth),]$density - 
        Subdata2[Subdata2$depth == min(Subdata2$depth),]$density)/
  mean(Subdata2$density)
st3<-(Subdata3[Subdata3$depth == max(Subdata3$depth),]$density - 
        Subdata3[Subdata3$depth == min(Subdata3$depth),]$density)/
  mean(Subdata3$density)
mean(st1)
mean(st2)
mean(st3)

######################
### Strat index 2 and mixed layer depth
##########################



require(plyr)

k=19#monthy
k2=3 #station
i=as.numeric(levels(as.factor(Datasum$month)))
Subdata<-Datasum[Datasum$station==k2 & Datasum$month == i[k],]
#Subdata$depth
Subdata$depth2<-round_any(Subdata$depth, 0.5)
SubdataAgg<-aggregate(. ~depth2, data=Subdata, mean)
j<-c(1:length(SubdataAgg$depth))
for (Val in j){
  SubdataAgg$DIFF<-SubdataAgg[j+1,6]-SubdataAgg[j,6]
}
#max(SubdataAgg$DIFF[3:max(j)-1])
SubdataAgg[SubdataAgg$DIFF == max(abs(SubdataAgg$DIFF[3:max(j)-1])),c(1,9)]
plot(SubdataAgg$depth ~ SubdataAgg$DIFF, ylim=c(40,0), xlim=c(-1,1))
abline(v=mean(SubdataAgg$DIFF[20:max(j)-1]))
abline(h=SubdataAgg[SubdataAgg$DIFF == max(abs(SubdataAgg$DIFF[3:max(j)-1])),9])
abline(h=1.5,
         col="red")
sort(SubdataAgg$DIFF)



######################
### integrate Chl over water column
##########################



Datasum<-read.table("All_CTD_Rampro_Chl.txt", header=TRUE)
Datasum$Fcal <- Datasum$F
Datasum$Fcal[Datasum$station == 1 & Datasum$month >1.5]<-Datasum$Fcal[Datasum$station == 1 & Datasum$month >1.5]*4.0585-0.2033
Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month != 7.47]<-Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month != 7.47]*4.0585-0.2033
Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month != 7.47]<-Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month != 7.47]*4.0585-0.2033
Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month == 7.47]<-Datasum$Fcal[Datasum$station == 2 & Datasum$month >1.5 & Datasum$month == 7.47]*0.2782+0.0461
Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month == 7.47]<-Datasum$Fcal[Datasum$station == 3 & Datasum$month >1.5 & Datasum$month == 7.47]*0.2782+0.0461
Datasum$Fcal[Datasum$month < 1.5]<-Datasum$Fcal[Datasum$month < 1.5]
str(Datasum)


k=11 #month
k2=1 #station
i=as.numeric(levels(as.factor(Datasum$month)))
Subdata<-Datasum[Datasum$station==k2 & Datasum$month == i[k],]
#Subdata$depth
Subdata$depth2<-round_any(Subdata$depth, 0.5)
SubdataAgg<-aggregate(. ~depth2, data=Subdata, mean)
j<-c(1:length(SubdataAgg$depth))
#for (Val in j){
# # SubdataAgg$Fcal<-SubdataAgg[j+1,12]-SubdataAgg[j,12]
#}
#max(SubdataAgg$DIFF[3:max(j)-1])
SubdataAgg$Fcal <- SubdataAgg$Fcal*0.5
SubdataAgg$Fcal[SubdataAgg$Fcal <0] <- 0
sum(SubdataAgg$Fcal)
SubdataAgg$depth[SubdataAgg$Fcal == max(SubdataAgg$Fcal)]

plot(SubdataAgg$depth ~ SubdataAgg$Fcal, ylim=c(40,0), type="b")
abline(h=SubdataAgg$depth[SubdataAgg$Fcal == max(SubdataAgg$Fcal)], col="green", lty=2)
abline(v=0)
abline(h=0)
#SubdataAgg$Fcal






#######

Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)

Dens<-Datasum[Datasum$station==3 & Datasum$month == 6.42,]$density
Dep<-Datasum[Datasum$station==3 & Datasum$month == 6.42,]$depth
Dens<-Dens[1:which(Dep==max(Dep))]
Dep<-Dep[1:which(Dep==max(Dep))]

f<-function(z) {(Pavg - P) * 9.81}

x <- Dep
P <- Dens
Pavg <- mean(Dens)

integrate(f, lower=-(max(Dep)), upper=0)

#######

#####################################################################
### Read out specific data

Datasum<-read.table("All_CTD_Rampro.txt", header=TRUE)


str(Datasum)
levels(as.factor(Datasum$month))
levels(as.factor(Datasum$station))
# Depth: 5m, 30m, 45m (st1), 75m (st2), 100m (st3)
Depth<-c(5,30,45,75,100)
Datasum[Datasum$month <2.7,]

Datasum<-Datasum[,c(1:3,6,7,9,10)]

aggregate(. ~station+month, data=Datasum[Datasum$depth > 4.8 & Datasum$depth <5.2, ], mean)
aggregate(. ~station+month, data=Datasum[Datasum$depth > 29.8 & Datasum$depth <30.2, ], mean, na.rm=TRUE)
aggregate(. ~station+month, data=Datasum[Datasum$depth > 44.5 & Datasum$depth <45.5 & Datasum$station == 1, ], mean, na.rm=TRUE)
aggregate(. ~station+month, data=Datasum[Datasum$depth > 74.8 & Datasum$depth <75.2 & Datasum$station == 2, ], mean, na.rm=TRUE)
aggregate(. ~station+month, data=Datasum[Datasum$depth > 99.8 & Datasum$depth <100.2 & Datasum$station == 3, ], mean, na.rm=TRUE)

max(Datasum[Datasum$month == 15.1 & Datasum$station == 1,]$depth)
aggregate(. ~station+month, data=Datasum[Datasum$depth > 34  & Datasum$station == 1, ], mean, na.rm=TRUE)

aggregate(. ~station+month, data=Datasum[Datasum$depth > 29 & Datasum$depth <31 &
                                          Datasum$station == 1, ], mean, na.rm=TRUE)


aggregate(. ~station+month, data=Datasum[Datasum$depth > 99.8 & Datasum$depth <100.2 & Datasum$station == 3, ], mean, na.rm=TRUE)

