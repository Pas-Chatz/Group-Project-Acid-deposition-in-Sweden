####Libraries
library(data.table)
library(ggplot2)
library(reshape2)
library(vegan)
###DAta
pH=read.csv(list.files(recursive = TRUE,pattern="pH",full.names = TRUE),
            fileEncoding ="UTF-8-BOM",stringsAsFactors = FALSE)
Phyto_1993=read.csv(list.files(recursive = TRUE,pattern="phytoplankton_1993",full.names = TRUE),
                    fileEncoding ="UTF-8-BOM",stringsAsFactors = FALSE)
Phyto_2010=read.csv(list.files(recursive = TRUE,pattern="phytoplankton_2010",full.names = TRUE),
                    fileEncoding ="UTF-8-BOM",stringsAsFactors = FALSE)

#Sampled lakes
TableOfSampleYears1993=table(Phyto_1993$Lake,Phyto_1993$Year)
TableOfSampleYears2010=table(Phyto_2010$Lake,Phyto_2010$Year)

FullSampling1993=TableOfSampleYears1993[rowSums(TableOfSampleYears1993==0)<=2,]#<=1 (keeps also sampling that miss one year)
FullSampling2010=TableOfSampleYears2010[rowSums(TableOfSampleYears2010==0)<=2,]

FullSampledLakes1993and2010=rownames(FullSampling1993)[rownames(FullSampling1993)%in%rownames(FullSampling2010)==TRUE]

FullSampledLakes1993and2010



####
lakes=split(Phyto_1993[Phyto_1993$Lake%in%FullSampledLakes1993and2010,],Phyto_1993$Lake[Phyto_1993$Lake%in%FullSampledLakes1993and2010])#only for 1993, split lake that are in common
lakes=lapply(lakes,function(q){
  aa=dcast(data = q,Year~taxon_name,value.var = "Biovolume",function(u){
    sum(u,na.rm=TRUE)
  },fill = 0)
  rownames(aa)=aa$Year
  aa=aa[,-1]
  aa=aa[rowSums(aa)!=0,]
  return(aa)
  
  
})
####CHARLIE ADD I addded this in. I was thinkning the other day that whilst 
ORDS=lapply(lakes,function(q){
  metaMDS(q)
  }
  )

pH=pH[pH$Lake%in%FullSampledLakes1993and2010,]

Ef=split(names(ORDS),names(ORDS))

for(i in 1:length(Ef)){
  nam=Ef[[i]]
  or=ORDS[[nam]]
  env=pH[pH$Lake==nam&pH$Year%in%as.numeric(rownames(or$points)),]
  env=aggregate(data=env,pH~Year,mean)
  if(nrow(env)>1){
  ef=envfit(or,env$pH)
  }else{
    NULL
  }
  Ef[[i]]=ef
}

summ_lakes=c()
for(i in 1:length(Ef)){
  nam=names(Ef)[i]
  ef=Ef[[i]]
  p=ef$vectors$pvals
  r2=ef$vectors$r
  axes=ef$vectors$arrows
attributes(axes)=NULL
summ_lakes=rbind(summ_lakes,  
data.frame(lake=nam,
           NMDS1=axes[1],
             NMDS2=axes[2],
             r2=r2,
             p=p,
             perm=ef$vectors$permutations))
  }
summ_lakes[summ_lakes$p<0.05,]

sapply(lakes,ncol) #richness per lake


as.matrix(dst$`Dagskärsgrund N`)[,1]
as.matrix(dst$Fagertärn)[,1]
as.matrix(dst$`Granfj. Djurgårds Udde`)[,1]

as.matrix(dst$`Megrundet N`)[,1]

as.matrix(dst$`S. Björkfjärden SO`)[,1]

as.matrix(dst$Spjutsjön)[,1]

as.matrix(dst$Stensjön)[,1]

as.matrix(dst$`Stora Envättern`)[,1]

as.matrix(dst$Tärnan)[,1]
as.matrix(dst$`Tärnan SSO`)[,1]





lakes2010=split(Phyto_2010[Phyto_2010$Lake%in%FullSampledLakes1993and2010,],Phyto_2010$Lake[Phyto_2010$Lake%in%FullSampledLakes1993and2010])#only for 2010, split lake that are in common
lakes2010=lapply(lakes2010,function(q){
  aa=dcast(data = q,Year~TaxonId,value.var = "Biovolume",function(u){
    sum(u,na.rm=TRUE)
  },fill = 0)
  rownames(aa)=aa$Year
  aa=aa[,-1]
  return(aa)
  
  
})
dst=lapply(lakes2010,vegdist)#add$ we can see every single lake
ord=metaMDS(lakes2010$`Ekoln Vreta Udd`)
ordiplot(ord,type="t",display = "site")    
