library(ggplot2)
library(dplyr)
library(data.table)
library(vegan)
library(reshape2)
##DAta
ph2010<-read.csv ("ph_2010-2018_UTF.csv", sep = ",",encoding = "UTF-8",stringsAsFactor=FALSE)
names(ph2010)[names(ph2010)=="X.U.FEFF.Lake"] <- "Lake"
ph1993<-read.csv ("pH_1993-2001_original_UTF8.csv", sep = ",",encoding = "UTF-8",stringsAsFactor=FALSE)
phyt1993<-read.csv("phytoplankton_1993_2001_data.csv", header = T, encoding = "UTF-8")
phyt2010<-read.csv("phytoplankton_2010_2018_data.csv", header = T, encoding = "UTF-8")
names(phyt2010)[names(phyt2010)=="X.U.FEFF.Lake"] <- "Lake"
names(phyt1993)[names(phyt1993)=="X.U.FEFF.Lake"] <- "Lake"
phyt2010$Biovolume <- as.numeric(as.character(phyt2010$Biovolume))
##########Sampled lakes for phytoplankton datase
(TableOfSampleYears1993=table(phyt1993$Lake,phyt1993$Year))
(TableOfSampleYears2010=table(phyt2010$Lake,phyt2010$Year))
(FullSampling1993=TableOfSampleYears1993[rowSums(TableOfSampleYears1993==0)<=2,])#<=1 (keeps also sampling that miss one year)
(FullSampling2010=TableOfSampleYears2010[rowSums(TableOfSampleYears2010==0)<=2,])
FullSampledLakes1993and2010=rownames(FullSampling1993)[rownames(FullSampling1993)%in%rownames(FullSampling2010)==TRUE]
FullSampledLakes1993and2010
##subset of which lakes are sampled in both time series
phyt1993_sampled<-phyt1993[phyt1993$Lake%in%FullSampledLakes1993and2010==TRUE,] #in those sampled years
phyt2010_sampled<-phyt2010[phyt2010$Lake%in%FullSampledLakes1993and2010==TRUE,]
###########for dataset2010
#i subset for month July-Aug
data_month_2010<-ph2010[ph2010$Month>=7&ph2010$Month<=8,]
data_month_1993<-ph1993[ph1993$Month>=7&ph1993$Month<=8,]
#subset for max depth <= 5m
data_depth_2010<-data_month_2010[data_month_2010$max_depth<=5,]
data_depth_1993<-data_month_1993[data_month_1993$max_depth<=5,]
####
#So far we have a subset of phytoplankton with sampled lakes 
#and for ph for Jul-Aug, max depth<=5 and continues year sampling. 
#Now i want to combine the lakes under all these conditions for all datasets.

##subset ph dataset lakes based on phytoplankton sampled lakes
TableOfSampleYears1993_ph=table(data_depth_1993$Lake,data_depth_1993$Year)
TableOfSampleYears1993_ph
TableOfSampleYears2010_ph=table(data_depth_2010$Lake,data_depth_2010$Year)
TableOfSampleYears2010_ph
FullSampling1993_ph=TableOfSampleYears1993_ph[rowSums(TableOfSampleYears1993_ph==0)<=2,]
FullSampling2010_ph=TableOfSampleYears2010_ph[rowSums(TableOfSampleYears2010_ph==0)<=2,]
FullSampledLakes1993and2010_ph=rownames(FullSampling1993_ph)[rownames(FullSampling1993_ph)%in%rownames(FullSampling2010_ph)==TRUE]
FullSampledLakes1993and2010_ph
####
ph2010<-data_depth_2010[data_depth_2010$Lake%in%FullSampledLakes1993and2010_ph==TRUE,] #in those sampled years
ph1993<-data_depth_1993[data_depth_1993$Lake%in%FullSampledLakes1993and2010_ph==TRUE,]
phyt1993_ph<-phyt1993_sampled[phyt1993_sampled$Lake%in%FullSampledLakes1993and2010_ph==TRUE,] #in those sampled years
phyt2010_ph<-phyt2010_sampled[phyt2010_sampled$Lake%in%FullSampledLakes1993and2010_ph==TRUE,]
#i return the previous argument to the phytoplankton dataset
ph1993_phyt<-ph1993[ph1993$Lake%in%FullSampledLakes1993and2010==TRUE,] #in those sampled years
ph2010_phyt<-ph2010[ph2010$Lake%in%FullSampledLakes1993and2010==TRUE,]
##### All include 87 lakes
#unique(ph1993_phyt$Lake)
#unique(ph2010_phyt$Lake)
#unique(phyt1993_ph$Lake)
#unique(phyt2010_ph$Lake)
unique(ph1993_phyt$Lake)%in%unique(ph2010_phyt$Lake) #this looks at the lakes present in 1993 to see if present in 2010
unique(ph2010_phyt$Lake)%in%unique(ph1993_phyt$Lake) #and vice versa
#Create mean pH
agg1993<-aggregate(data=ph1993_phyt, pH~Lake+Year, function(q)
{c(mean = mean(q), SD = sd(q))})
agg1993
agg1993<-cbind(agg1993[,c("Lake","Year")],agg1993$pH)
names(agg1993)[names(agg1993)=="mean"] <- "ph"
sapply(agg1993,class)
#take the mean of pH of months by lake.
agg2010<-aggregate(data=ph2010_phyt, pH~Lake+Year, function(q)
{c(mean = mean(q), SD = sd(q))})
agg2010<-cbind(agg2010[,c("Lake","Year")],agg2010$pH)
sapply(agg2010, class) # i cannot understand whats happening with the class of pH.
names(agg2010)[names(agg2010)=="mean"] <- "ph"
#combined dataset of ph values!!!
data.ph.sub = rbind(agg1993, agg2010)#combines two identical data frames 
names(data.ph.sub)
head(data.ph.sub)
##########################
#1993
abundance1993<-phyt1993_ph[c(1,19,30,36:37)]
abundance2010<-phyt2010_ph[c(1,19,30,36:37)]
abundance2010$Biovolume

summary(abundance1993$Biovolume) #i inspect the rance of biovolume to know how to subset
abundance1993$level<-NA
abundance1993$level<-"Rare"
abundance1993$level[abundance1993$Biovolume>=10]<-"Common"
#2010
summary(abundance2010$Biovolume)
abundance2010$level<-NA
abundance2010$level<-"Rare"
abundance2010$level[abundance2010$Biovolume>=1&abundance2010$Biovolume<10]<-"intermidiate"
abundance2010$level[abundance2010$Biovolume>=10]<-"Common"
#
ggplot(abundance1993, aes(fill=group,y=Biovolume, x=reorder(Biovolume,taxon_name))) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_blank())+
  ggtitle("Rank abundance of phytoplankton 1993") +
  xlab("Species") + ylab("Abundance (Biovolume)")+
  facet_wrap(~level, scale="free")

#SAME graph with r base introduced by Bob Muscarella
#1993
rank.abund_1993 <- rev(sort(tapply(abundance1993$Biovolume, factor(abundance1993$taxon_name), sum)))
barplot(rank.abund_1993, las=2, ylab="Biovolume", main="Rank Abundance of phytoplankton 1993")
#2010 something is wrong, still doesnt treat biovolume for 2010 as numeric
rank.abund_2010<- rev(sort(tapply(abundance2010$Biovolume,
                                  factor(abundance2010$taxon_name), sum))) 
barplot(rank.abund_2010, las=2, ylab="Biovolume", main="Rank Abundance of phytoplankton 2010")
###Compute species richness
#Compute species richness for phyt (i.e., number of taxon_name per Lake)_Bob  
sp.rich_1993 <- rowSums(table(phyt1993_ph$Lake, phyt1993_ph$taxon_name)>0)
sp.rich_1993
sp.rich_1993 <- data.frame(Lake=names(sp.rich_1993), richness = sp.rich_1993)
sp.rich_2010<-rowSums(table(phyt2010_ph$Lake, phyt2010_ph$taxon_name)>0)

sp.rich_2010 <- data.frame(Lake=names(sp.rich_2010),richness = sp.rich_2010)

sp.rich_total<-merge(sp.rich_1993,sp.rich_2010, by.x = "Lake", by.y = "Lake")
unique(sp.rich_total$Lake) 
#####
names(sp.rich_total)[names(sp.rich_total)=="richness.x"] <- "1993-2001"
names(sp.rich_total)[names(sp.rich_total)=="richness.y"] <- "2010-2018"
sp.rich_total$dif<-(sp.rich_total$`2010-2018` - sp.rich_total$`1993-2001`)
#
sp.rich_total$change_prer<-(sp.rich_total$dif/sp.rich_total$`1993-2001`)*10
sp.rich_total
sp.rich_total$Change <- NA
sp.rich_total$Change <- "Decreased or constant"
sp.rich_total$Change[sp.rich_total$dif > 0 & sp.rich_total$dif <= 20] <- "Increased"
sp.rich_total$Change[sp.rich_total$dif > 20] <- "Increased a lot"
sp.rich_total
#
#i transform to long format in order to use ggplot
sp.rich.melted <- melt(sp.rich_total,
                       id=c("Lake","dif","Change","change_prer"),
                       variable.name = "Year",
                       value.name = "Richness")
#i change the order to alphabetic
sp.rich.melted<-sp.rich.melted[order(sp.rich.melted$Lake),]
sp.rich.melted
# 
ggplot(sp.rich.melted, aes(fill=Year, y=Richness, x=reorder(Lake,Richness))) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.text.x=element_blank())+
  ggtitle("Species richness per lake") +
  xlab("Lake") + ylab("Species Richness") + facet_wrap(~Change, scale="free")

######### INSERTED BY BOB #########
#AND TRANSFORMED BY PASCHALIS
# Make vector of all the different lakes during both times
#i want to produce a data frame with species richness for each lake each year in order to correlate it
#with ph and aluminium for each lake each year!!!!
all_lakes <- unique(c(as.character(sp.rich_2010$Lake), 
                      as.character(sp.rich_1993$Lake)))
all_lakes
#WHY SO MANY LAKES? (Paschalis)
# Make a vector for years (one each of 1993 and 2010 for each lake)
year_1993 <- rep(c(1993,2001), times=length(all_lakes))
year_2010 <- rep(c(2010,2018), times=length(all_lakes)) #paschalis=i dont know how to proceed!
# Match richness (from above) with lake names
rich_1993 <- sp.rich_1993$richness[match(all_lakes, sp.rich_1993$Lake)]
rich_2010 <- sp.rich_2010$richness[match(all_lakes, sp.rich_2010$Lake)]

# Interleave richness
richness_1993_2010 <- c(rbind(rich_1993, rich_2010))
richness_1993_2010
# Match and interleave pH values
pH_1993_2010 <- c(rbind(agg1993$ph[match(all_lakes, ph1993_phyt$Lake)], 
                        agg2010$ph[match(all_lakes, ph2010_phyt$Lake)]))

mydf <- data.frame(lake = rep(all_lakes, each=2),
                   year = year,
                   richness = richness_1993_2010,
                   ph = pH_1993_2010)
#################################



