library(cowplot)
library(data.table)
library(dplyr)
library(colorRamps)
library(ggplot2)
library(reshape2)
library(vegan)
library(heatmap3)
library(ggpubr)
theme_set(theme_cowplot())

####Ph Data#### loads and subsets pH data and gets one data frame for pH
#2010-2018 dataset
ph2010<-read.csv ("ph_2010-2018_UTF.csv", sep = ",",fileEncoding =  "UTF-8-BOM",stringsAsFactor=FALSE)
#for some reason lake variable has to be changed
names(ph2010)[names(ph2010)=="X.U.FEFF.Lake"] <- "Lake"
names(ph2010)[which(grepl("Stationskoordinat.E.Y|Stationskoordinat.N.X",names(ph2010)))]=c("Y","X")
#1993-2001 dataset
ph1993<-read.csv ("pH_1993-2001_original_UTF8.csv", sep = ",",fileEncoding  = "UTF-8-BOM",stringsAsFactor=FALSE)
names(ph1993)[which(grepl("Stationskoordinat.E.Y|Stationskoordinat.N.X",names(ph1993)))]=c("Y","X")

ph2010<-ph2010[ph2010$Month>=7&ph2010$Month<=8&ph2010$max_depth<=5,]


##same story for 1993 dataset
#i subset for month June-July-Aug
ph1993<-ph1993[ph1993$Month>=7&ph1993$Month<=8&ph1993$max_depth<=5,]

TableOfSampleYears1993_ph=table(ph2010$Lake,ph2010$Year)
TableOfSampleYears1993_ph
TableOfSampleYears2010_ph=table(ph1993$Lake,ph1993$Year)
TableOfSampleYears2010_ph
FullSampling1993_ph=TableOfSampleYears1993_ph[rowSums(TableOfSampleYears1993_ph==0)<=2,]
FullSampling2010_ph=TableOfSampleYears2010_ph[rowSums(TableOfSampleYears2010_ph==0)<=2,]
FullSampledLakes1993and2010_ph=rownames(FullSampling1993_ph)[rownames(FullSampling1993_ph)%in%rownames(FullSampling2010_ph)==TRUE]
FullSampledLakes1993and2010_ph

ph2010<-ph2010[ph2010$Lake%in%FullSampledLakes1993and2010_ph==TRUE,] #in those sampled years
ph1993<-ph1993[ph1993$Lake%in%FullSampledLakes1993and2010_ph==TRUE,]

####Aluminium sampling  Same as pH for AL

AL93=read.csv("Al_1993-2001.csv",sep=",",fileEncoding = "UTF-8-BOM",stringsAsFactor=FALSE)
AL93=AL93[AL93$Month%in%c(7,8)&AL93$max_depth<=5,]
AL10=read.csv("Al_2010-2018.csv",sep=",",fileEncoding = "UTF-8-BOM",stringsAsFactor=FALSE)
AL10=AL10[AL10$Month%in%c(7,8)&AL10$max_depth<=5,]

ALsamples1993=table(AL93$Lake,AL93$Year)
ALsamples2010=table(AL10$Lake,AL10$Year)
ALsamples1993=ALsamples1993[rowSums(ALsamples1993==0)<=2,]
ALsamples2010=ALsamples2010[rowSums(ALsamples2010==0)<=2,]

ALLakes=rownames(ALsamples1993)[rownames(ALsamples1993)%in%rownames(ALsamples2010)]

####pH and Aluminium lakes#### 
ALLakes=ALLakes[ALLakes%in%FullSampledLakes1993and2010_ph]


####PhytoData#### Loads data for lake with only 2 years missing
phyt1993<-read.csv("phytoplankton_1993_2001_data.csv", header = T, fileEncoding =  "UTF-8-BOM")
phyt2010<-read.csv("phytoplankton_2010_2018_data.csv", header = T, fileEncoding =  "UTF-8-BOM")

phyt1993<-phyt1993[phyt1993$Month>=7&phyt1993$Month<=8&phyt1993$max_depth<=5,]
phyt2010<-phyt2010[phyt2010$Month>=7&phyt2010$Month<=8&phyt2010$max_depth<=5,]

phytsamplng93=rowSums(table(phyt1993$Lake,phyt1993$Year)==0)
phytsamplng93=names(phytsamplng93[phytsamplng93<=2])
phytsamplng10=rowSums(table(phyt2010$Lake,phyt2010$Year)==0)
phytsamplng10=names(phytsamplng10[phytsamplng10<=2])

Phyt9310=phytsamplng10[phytsamplng10%in%phytsamplng93]

LAKES=FullSampledLakes1993and2010_ph[FullSampledLakes1993and2010_ph%in%Phyt9310] #list of lakes

phyt1993=phyt1993[phyt1993$Lake%in%LAKES,]
phyt2010=phyt2010[phyt2010$Lake%in%LAKES,]

#converts <0.00005 to lowest value
phyt2010$Biovolume=as.numeric(as.character(phyt2010$Biovolume))       
phyt2010=phyt2010[is.na(phyt2010$Biovolume)==FALSE,]

##creates list of separate decades with species names Genres removed grep selects only those vector positions that comply to the first arguement
phyt=list(phyt93=phyt1993[grep(" ",phyt1993$taxon_name),],
          phyt10=phyt2010[grep(" ",phyt2010$taxon_name),])

###Aggregate across years per decade
phyt_agg=lapply(phyt,function(q){
  aggregate(data=q,Biovolume~Lake+taxon_name+TaxonId,mean)
})


##community table across decades 
phyt_comm=lapply(phyt_agg,function(q){
  aa=dcast(data = q,Lake~taxon_name,value.var = "Biovolume",function(u){
    sum(u,na.rm=TRUE)
  },fill = 0)
  row.names(aa)=aa$Lake
  aa=aa[,-1]
  return(aa)
})

ords=lapply(phyt_comm,metaMDS)
heatmap3::heatmap3(t(phyt_comm[[1]]),distfun = function(q){vegdist(q)},
                   method="average")

heatmap3::heatmap3(t(phyt_comm[[2]]),distfun = function(q){vegdist(q)},
                   method="average")

spMissing=colnames(phyt_comm$phyt93)[colnames(phyt_comm$phyt93)%in%colnames(phyt_comm$phyt10)==FALSE]
spAppearing=colnames(phyt_comm$phyt10)[colnames(phyt_comm$phyt10)%in%colnames(phyt_comm$phyt93)==FALSE]
spRemaining=colnames(phyt_comm$phyt93)[colnames(phyt_comm$phyt93)%in%colnames(phyt_comm$phyt10)==TRUE]



##species richness, abundace and diversityb measures across decades
sp_rich=lapply(phyt_comm, function(q){
  data.frame(sp_richness=rowSums(q>0),
             biovolme=rowSums(q),
             shannon=diversity(q,index="shannon"),
             simpson=diversity(q,"simpson"))
})

#bind list together
sp_richness=rbind(data.frame(sp_rich$phyt93,dec=factor('1993 - 2001')),
                  data.frame(sp_rich$phyt10,dec=factor('2010 - 2018')))

##plots NOTE SOME ARE ON A LOG SCALE
ggplot(data=sp_richness,aes(x=dec,y=sp_richness, fill=dec))+
  geom_boxplot()+
  labs(title="Species richness comparison between time series",x="Time Series", y = "Species Richness")+
  labs(fill = "Time Series")

# means and medians for each time series of sp_richness 
sp_richness %>%
  group_by(dec) %>%
  summarise(sp_richness = mean(sp_richness)) -> mean_sprich
mean_sprich
sp_richness %>%
  group_by(dec) %>%
  summarise(sp_richness = median(sp_richness)) -> summ_sprich
summ_sprich
##calculate the mean sp.richness per index per decade. (Paschalis)
mean(sp_richness$shannon[sp_richness$dec=='1993 - 2001'])
mean(sp_richness$shannon[sp_richness$dec=='2010 - 2018'])
mean(sp_richness$simpson[sp_richness$dec=='1993 - 2001'])
mean(sp_richness$simpson[sp_richness$dec=='2010 - 2018'])
## not evenly distributed 

ggplot(data=sp_richness,aes(x=dec,y=biovolme, fill=dec))+
  geom_boxplot()+scale_y_log10()+
  labs(title="Biovolume comparison between time series",x="Time Series", y = "log(Biovolume)")+
  labs(fill = "Time Series")

mod=lm(data=sp_richness,sp_richness~dec)
mod=lm(data=sp_richness,biovolme~dec)

anova(mod)
ggplot(data=sp_richness,aes(x=sp_richness,y=biovolme,color=dec))+
  geom_point()+
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,fullrange=TRUE)+
  #facet_wrap(~dec,scale="free")+
  scale_y_log10()+
  ylab("Biovolume (UNITS)")+
  xlab("Species Richness")+
  labs(colour = "Time Series")

sha<-ggplot(data=sp_richness,aes(x=dec,y=shannon, fill=dec))+
  geom_boxplot()+
  labs(title="Shannon Diversity",x="Time Series", y = "Shannon's index")

simp<-ggplot(data=sp_richness,aes(x=dec,y=simpson, fill=dec))+
  geom_boxplot()+
  labs(title="Simpson Diversity",x="Time Series", y = "Simpson's Diversity")+
  labs(fill = "Time Series")
simp + theme(legend.position="none")
sha + theme(legend.position="none")
ggarrange(simp,sha, ncol=2, nrow=1, common.legend = TRUE, legend="right")

grid.arrange(sha, simp, nrow = 1)


###rank abundance curves per lake. It may be worth looking to see what is at number 1 in seprarate decades
RankAbund=lapply(phyt,function(q){
  xx=lapply(LAKES,function(x){
    aa=dcast(data = q[q$Lake==x,],Lake~taxon_name,value.var = "Biovolume",function(u){
    sum(u,na.rm=TRUE)
  },fill = 0)
    dat=data.frame(lake=x,
                   taxon=colnames(aa[,-1]),
                   abund=colSums(aa[,-1]),
                   rank=rank(-colSums(aa[,-1])))
    return(dat)
    
    
})
dat=do.call("rbind",xx)
return(dat)
})
##curve plots
RankAbund=rbind(data.frame(RankAbund[[1]],dec=factor(1993)),
                data.frame(RankAbund[[2]],dec=factor(2010)))
#i change the names (Paschalis)
levels(RankAbund$dec)[levels(RankAbund$dec)=='1993'] <- "1993-2001"
levels(RankAbund$dec)[levels(RankAbund$dec)=='2010'] <- "2010-2018"

ggplot(data=RankAbund,aes(x=rank,y=abund,colour=lake))+
  geom_point()+
  facet_wrap(~dec)+
  scale_y_log10()+
  theme(legend.position = "none")+
  ylab("Abundance")+
  xlab("Rank")


####Beta diversity between decades
dist=c()
for(i in seq(length(LAKES))){
  LAK=LAKES[i]
  com=merge(phyt_comm$phyt93[rownames(phyt_comm$phyt93)==LAK,],
  phyt_comm$phyt10[rownames(phyt_comm$phyt10)==LAK,],all.x=TRUE,all.y=TRUE)
  com[is.na(com)]=0
  
  B=vegdist(com)
  PA=com
  PA[PA>0]=1
  
  PAd=vegdist(PA)
  beta=betadiver(com,"w")
  dist=rbind(dist,data.frame(Lake=LAK,BCdist=B[1],PAdist=PAd[1],beta=beta[1]))
}

####pH differences#### finds differences between lakes in pH can be converted to AL


phagg2010=aggregate(data=ph2010[ph2010$Lake%in%LAKES,],pH~Lake,function(q){pH2010=mean(q)})
phagg1993=aggregate(data=ph1993[ph1993$Lake%in%LAKES,],pH~Lake,function(q){pH1993=mean(q)})
names(phagg2010)=c("Lake","pH2010")
names(phagg1993)=c("Lake","pH1993")

phagg=merge(phagg1993,phagg2010)

phagg$pHdiff=phagg$pH2010-phagg$pH1993

phagg$Lake==rownames(sp_rich$phyt93)
phagg$Lake==rownames(sp_rich$phyt10)

phagg$SPr93=sp_rich$phyt93$sp_richness
phagg$SPr10=sp_rich$phyt10$sp_richness

phagg$SPrdiff=phagg$SPr10-phagg$SPr93

phagg=merge(phagg,dist)

ggplot(data=phagg[phagg$pHdiff!=max(phagg$pHdiff),],aes(x=pHdiff,y=beta))+  #the !=max removes the max outlier of pH
  geom_point()


summary(mod)
ggplot(data=phagg,aes(x=SPrdiff,y=SPrdiff))+
  geom_point()

locations=merge(phagg,ph1993[,c(1,4:5)])
locations=unique(locations)
head(locations)
ggplot(data=locations,aes(x=X,y=Y,colour=beta))+
  geom_point()+
  scale_color_gradientn(colours = matlab.like2(100))

#QGIS
coord<-write.csv(locations,"coord.csv", row.names = FALSE)
