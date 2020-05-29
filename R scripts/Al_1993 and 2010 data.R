library(ggplot2)
library(dplyr)
library(data.table)
library(vegan)
library(reshape2)
Al2010<-read.csv("Al_2010-2018.csv", sep = ",",encoding = "UTF-8")
Al1993 <- read.csv("Al_1993-2001.csv", sep = ",", fileEncoding = "UTF-8-BOM")
str(Al2010)
names(Al2010)[names(Al2010)=="X.U.FEFF.Lake"] <- "Lake"
names(Al1993)[names(Al1993)=="X.U.FEFF.Lake"] <- "Lake"

TableOfSampleYears1993_Al=table(Al2010$Lake,Al2010$Year)
TableOfSampleYears1993_Al
TableOfSampleYears2010_Al=table(Al1993$Lake,Al1993$Year)
TableOfSampleYears2010_Al
FullSampling1993_Al=TableOfSampleYears1993_Al[rowSums(TableOfSampleYears1993_Al==0)<=2,]
FullSampling2010_Al=TableOfSampleYears2010_Al[rowSums(TableOfSampleYears2010_Al==0)<=2,]
FullSampledLakes1993and2010_Al=rownames(FullSampling1993_Al)[rownames(FullSampling1993_Al)%in%rownames(FullSampling2010_Al)==TRUE]
FullSampledLakes1993and2010_Al

Al2010<-Al2010[Al2010$Lake%in%FullSampledLakes1993and2010_Al==TRUE,] #in those sampled years
Al1993<-Al1993[Al1993$Lake%in%FullSampledLakes1993and2010_Al==TRUE,]

#2010
#subset what variables you want
data2010<-Al2010[,c(1,4:8,18:23,26:31)]
sapply(data2010,class)
head(data2010)
data2010[8:18]<-lapply(data2010[8:18], as.numeric)
#i changed year, month, dat, mindepth, maxdepth to numeric
names(data2010)
sapply(data2010,class)

# subset for July and Aug
Aldata_month_2010<-data2010[data2010$Month>=7&data2010$Month<=8,]
#subset for max depth <= 5m
Aldata_depth_2010<-Aldata_month_2010[Aldata_month_2010$max_depth<=5,]
#continues years
Aldata.continuess_2010<-Aldata_depth_2010[Aldata_depth_2010$Year>=2010&Aldata_depth_2010$Year<=2018,]
all.equal(Aldata_depth_2010,Aldata.continuess_2010)
#Aldata.continuess_2010 contains the lakes whch have continues mesurements from 2010 to 2018.

sum(is.na(Aldata_depth_2010$Al_NAD)) 
# tells me how many NA values there are
# 480 for AL_ICPAES
# 957 for AL_NA
# 781 for Al_NAD



#take the mean of Al of months by lake.
Alagg2010<-aggregate(data=Aldata.continuess_2010, Al_NAD~Lake+Year, function(q)
{c(mean = mean(q), SD = sd(q))})
Alagg2010<-cbind(Alagg2010[,c("Lake","Year")],Alagg2010$Al_NAD)
sapply(Alagg2010, class) 
###



##same story for 1993 dataset
names(Al1993)
data1993<-Al1993[,c(1,4:8,18:23,26:30)]
sapply(data1993,class)
head(data1993)
names(data1993)
data1993[7:17]<-lapply(data2010[7:17], as.numeric)
#i changed year, month, dat, mindepth, maxdepth to numeric
names(data1993)
sapply(data1993,class)
#i subset for month June-July-Aug
Aldata_month_1993<-data1993[data1993$Month>=7&data1993$Month<=8,]
names(Aldata_month_1993)
#subset for max depth <= 5m
Aldata_depth_1993<-Aldata_month_1993[Aldata_month_1993$max_depth<=5,]
#continues years
Aldata.continuess_1993<-Aldata_depth_1993[Aldata_depth_1993$Year>=1993&Aldata_depth_1993$Year<=2001,]
all.equal(Aldata_depth_1993,Aldata.continuess_1993)
#take the mean of pH of moths by lake.
sum(is.na(Aldata.continuess_1993$Al_NAD)) 
# tells me how many NA values there are
## Al_NAD = 0
## no continuous observations hmmmmm

Alagg1993<-aggregate(data=Aldata.continuess_1993, Al_NAD~Lake+Year, function(q)
{c(mean = mean(q), SD = sd(q))})
Alagg1993<-cbind(Alagg1993[,c("Lake","Year")],Alagg1993$Al_NAD)

#combined dataset of ph values
Aldata.ph.sub = rbind(Alagg1993, Alagg2010)#combines two identical data frames 
  # only up to 2013
