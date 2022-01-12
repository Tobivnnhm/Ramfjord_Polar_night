library(readxl)
library(vegan)

setwd("D:/UiT/1UiT/artikel/r")

data <- as.data.frame(read_xlsx("R_PN_Barplot_GatheredData.xlsx"))
str(data)
row.names(data)<-data[,1]
data<-data[,-1]
data<-as.data.frame(t(data))
data$Month<-substr(row.names(data),start=1,stop=3)
data$Type<-substr(row.names(data),start=5,stop=7)
data$Station<-substr(row.names(data),start=nchar(row.names(data)), stop=nchar(row.names(data)))



#1) if the community structure is different between the months (calculating away the effects of sampling type)

adonis(data[,1:17] ~ Month + Type + Station, data=data, strata=data$Type, method="bray", permutations=999)



#2) if the sampling types are different (calculating away the effects of time).

adonis(data[,1:17] ~ Month + Type + Station + Month*Type, strata=data$Month, data=data, method="bray", permutations=999)

