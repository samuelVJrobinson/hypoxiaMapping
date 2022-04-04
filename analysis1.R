#FIRST APPROACH TO ANALYSIS: FIT MODELS OF DIFFERENT SATELLITE CHANNELS (SMOOTH) AND MODEL BOTTOM DO BASED ON THESE
#WRITTEN BY SR, SPRING 2021

# Load everything ---------------------------------------------------------
library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(ggmap)
library(ggpubr)
library(animation)
library(mgcv)
library(vegan)
library(gstat)
library(beepr)

setwd("~/Documents/hypoxiaMapping")

source('helperFunctions.R')

# Load data --------------------

#Note:
# sample_waterData.csv was used for analyses sets 1 and 2 (spring-fall 2021), sample_waterData3.csv used during spring 2022

#PROBLEM: the two water datasets are different from each other. Ask Yingjie

#Water data from scratch
wDat <- read.csv('./data/sample_waterData3.csv') %>% mutate(Date=as.Date(Date,'%Y-%m-%d')) %>%
  mutate(doy=as.numeric(format(Date,format='%j'))) %>%
  # mutate(YEID=ifelse(YEID=='2014_149','2014_003',YEID)) %>% #Set 2014_149 as 2014_003 (see below)
  group_by(YEID) %>% mutate(maxDepth=max(Depth,Depth_dem)) %>% mutate(propDepth=Depth/maxDepth) %>%
  select(-Depth_dem) %>% relocate(YEID:Date,doy,Depth,maxDepth,propDepth,DO,Temp:lon) %>% ungroup() %>%
  st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) #Louisiana offshore

#Check new data - a bunch in new locations
# wDat %>% mutate(isNew=grepl('hyp_watch',YEID)) %>% ggplot()+geom_sf(aes(col=isNew,size=isNew))

Emean <- mean(unique(wDat$E)); Nmean <- mean(unique(wDat$N)) #Center values of E/N, used for scaling

wDat <- wDat %>%
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  st_transform(4326)

# #All YEIDs are from a single day
# wDat %>% st_drop_geometry() %>% group_by(YEID) %>% summarize(nDays=length(unique(Date))) %>% filter(nDays!=1)

# # 2014_003, 2014_149 are different times at the same location
# wDat %>% mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2]) %>% unite(loc,x,y) %>% st_drop_geometry() %>% select(YEID,loc) %>%
#   distinct() %>% group_by(loc) %>% mutate(n=n()) %>% filter(n>1)

#Get locations from wDat to join onto spatial data
locIndex <- wDat %>% select(YEID,Date,doy,E:sN,geometry) %>% distinct()

#Spectral data
# path <- './data/sample_spectralData.csv' #Older dataset (March 2021)
# path <- './data/sample_spectralData2.csv' #Older dataset (September 2021)
path <- './data/sample_spectralData3.csv' #Newer dataset (March 2022)

sDat <- read.csv(path) 

#Check new values
sDat$YEID[!sDat$YEID %in% wDat$YEID] %>% unique #Entries in sDat not in wDat - OK
wDat$YEID[!wDat$YEID %in% sDat$YEID] %>% unique #OK
locIndex$YEID[!locIndex$YEID %in% sDat$YEID] %>% unique #locIndex entries not in sDat - OK
locIndex$YEID[!locIndex$YEID %in% wDat$YEID] %>% unique #locIndex entries not in wDat - OK

sDat <- sDat %>% 
  left_join(locIndex,by='YEID') %>% #Join in spatial info
  select(-Date,-doy) %>%
  st_as_sf() %>%
  mutate(date_img=as.Date(date_img,'%Y-%m-%d'),doy=as.numeric(format(date_img,'%j'))) %>%
  mutate(nflh=ifelse(nflh>(1-1e-5),(1-1e-5),nflh)) #Set upper limit of nflh just below 1

#Set lowest positive number for spectral variables
# lwrLimits <- c(0.1026964,1e-05,32.8,1e-09,1e-09,1.8001e-05,0.000452001,0.000770001,0.000822001,0.000566001,1e-09,1e-09,1e-09) #Older values
lwrLimits <- apply(st_drop_geometry(sDat)[,3:16],2,function(x) min(x[x>0],na.rm=TRUE)) #Minimum nonzero values
names(lwrLimits) <- c('chlor_a','nflh','poc','Rrs_412','Rrs_443','Rrs_469','Rrs_488','Rrs_531','Rrs_547','Rrs_555','Rrs_645','Rrs_667','Rrs_678','sst')

sDat <- sDat %>%
  #Rescales negative values be above 0.95*minimum positive value
  mutate(across(chlor_a:sst,~ifelse(.x<lwrLimits[names(lwrLimits)==cur_column()],lwrLimits[names(lwrLimits)==cur_column()]*0.95,.x))) %>%
  mutate(numNA=apply(st_drop_geometry(.)[,3:16],1,function(x) sum(is.na(x)))) %>% #Count NAs in data columns
  mutate(propNA=numNA/max(numNA))

#choose only surface water, using order of depth measurements
surfWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth)) %>% ungroup() %>%
  filter(depthOrd==1,Depth<80) %>% select(-depthOrd)

#choose only bottom water
bottomWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth,decreasing=TRUE)) %>% ungroup() %>%
  filter(depthOrd==1,Depth<80) %>% select(-depthOrd)

#PCA imputation of missing data imputation for missing data
# Examples: https://link.springer.com/article/10.1007/s11258-014-0406-z
# http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html

#Most measurements are missing 13 or 14 measurements
table(sDat$numNA)

#"Most missing" values is chlor_a. Rrs_ values are all missing together
apply(st_drop_geometry(select(sDat,chlor_a:sst)),2,function(x) sum(is.na(x)))/nrow(sDat)

#Looks like the missing values split into two clusters
hist(sDat$propNA,xlab='Proportion NA values')

#Get all values into a dataframe
sDatMat <- sDat %>% st_drop_geometry() %>%
  select(YEID:sst,propNA) %>%
  unite(ID,YEID,date_img) %>% remove_rownames() %>%
  column_to_rownames(var='ID') %>%
  filter(propNA<0.3) %>% #Drop rows that have more than 30% of their values missing
  select(-propNA) %>%
  mutate(across(everything(),log)) #Log-scale variables

# #For complete data (no NAs)
# sDatMat_noNA <- as.matrix(sDatMat[apply(sDatMat,1,function(x) sum(is.na(x)))==0,]) #No NA values
# pca <- prcomp(sDatMat_noNA,scale. = TRUE)
# plot(1:14,cumsum(pca$sdev^2)/sum(pca$sdev^2),xlab='PC',ylab='% Var',pch=19,type='b')
# abline(h=0.95,col='red') #Looks like about 6 dims are needed for 95% of variance
# #Singular value decomposition - identical
# svdMat <- svd(scale(sDatMat_noNA)) #First 2 dimensions
# plot(svdMat$d^2/sum(svdMat$d^2),ylab='Relative variance',xlab='Eigenvalue',pch=19)
# plot(cumsum(svdMat$d^2)/sum(svdMat$d^2),ylab='Cumulative variance',xlab='Eigenvalue',pch=19,type='b')
# abline(h=0.95,col='red') #Looks like about 6 dims are needed for 95% of variance
# plot(log(pca$sdev),log(svdMat$d)) #Similar things

library(missMDA)
# n_components <- estim_ncpPCA(sDatMat, verbose = TRUE,method="Regularized",method.cv="gcv",ncp.min=1,ncp.max=13) #Takes about 5 mins
# plot(1:13,n_components$criterion,xlab='Number of Dimensions',ylab='GCV Criterion',type='b',pch=19) #Looks like about 7 components is OK for prediction
sDatMat_imputed <- imputePCA(sDatMat,ncp=7,scale=TRUE,method='Regularized') #Impute missing data using 7 dimensions

head(sDatMat) #Original data
image(as.matrix(sDatMat))
head(sDatMat_imputed$completeObs) #Filled-in data
image(as.matrix(sDatMat_imputed$completeObs))

pca1 <- prcomp(sDatMat_imputed$completeObs,scale=TRUE) #Get PCA values

# Looks like about 6 dims are needed for 95% of variance
(p <- data.frame(pc=1:length(pca1$sdev),cVar=cumsum(pca1$sdev^2)/sum(pca1$sdev^2)) %>%
  ggplot(aes(x=pc,y=cVar))+geom_point()+geom_line()+
  labs(x='Principle Component',y='Cumulative Variance')+
  geom_hline(yintercept = 0.95,col='red',linetype='dashed'))
cumsum(pca1$sdev^2)/sum(pca1$sdev^2) #95% of var in first 6 PCs
ggsave('./figures/pcVar.png',p,width=5,height=5)

#Factor loadings of first 6 PCs
(p <- pca1$rotation[,1:6] %>% data.frame() %>% rownames_to_column(var='var') %>%
  pivot_longer(PC1:PC6) %>%
  mutate(name=factor(name,labels=paste0('PC',1:6,': ',round(pca1$sdev[1:6]^2/sum(pca1$sdev^2),3)*100,'% Variance'))) %>%
  mutate(var=factor(var,levels=rev(unique(var)))) %>%
  ggplot()+
  # geom_point(aes(y=var,x=value))+
  geom_col(aes(y=var,x=value))+geom_vline(xintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~name)+
  labs(x='Loading',y=NULL,title='Factor loadings for Principle Components 1-6'))
ggsave('./figures/factorLoadings.png',p,width=10,height=5)

#Make PC scores using center, scale, and rotation matrix
sDatMat_imputed$completeObs[1,] #Input data
pca1$x[1,] #PC scores (goal)
pca1$sdev #sqrt-Eigenvalues
pca1$rotation #Factor loadings
pca1$center #Mean of original data
pca1$scale #SD of original data

#Check the math...
((sDatMat_imputed$completeObs[1,]-pca1$center)/pca1$scale) %*% pca1$rotation #Reconstructed PCs
pca1$x[1,] #PCs from method
(((sDatMat_imputed$completeObs[1,]-pca1$center)/pca1$scale) %*% pca1$rotation)-pca1$x[1,] #Identical

#Join imputed PCA values back to original dataset (NA if no data on that day)
sDat_pca <- pca1$x[,1:6] %>% as.data.frame() %>% rownames_to_column('ID')
sDat <- sDat %>% unite(ID,YEID,date_img,sep='_',remove=FALSE) %>%
  left_join(sDat_pca,by='ID') %>% select(-ID) %>% mutate(gap=is.na(PC1))
rm(sDat_pca,sDatMat_imputed,sDatMat,p) #cleanup
save(bottomWDat,Emean,Nmean,locIndex,pca1,sDat,surfWDat,wDat,lwrLimits,file='./data/all2014_3.Rdata')

#Load data from saved file
load('./data/all2014_3.Rdata')

# Summary statistics of dataset -----------------------------

nrow(bottomWDat) #204 bottom DO measurements
length(unique(bottomWDat$doy)) #25 unique days
range(bottomWDat$Date) #Taken from May 26 - August 2
length(unique(with(bottomWDat,paste0(E,'_',N)))) #203 unique locations

# Take a look at water data ---------------------

#Base map
spDomain <- wDat %>% st_union() %>% st_convex_hull() %>% st_buffer(dist=0.25) %>% st_bbox() #Spatial domain
names(spDomain) <- c('left','bottom','right','top')
basemap <- get_map(location=spDomain) #Get from Google Earth

#Removing extra channels for now, and stripping out data-less days
sDat2 <- sDat %>% select(YEID,date_img,doy,chlor_a:poc,E:geometry) %>%
  mutate(noData=is.na(chlor_a)&is.na(nflh)&is.na(poc)) %>% #Strip days where all data are missing
  group_by(date_img) %>% mutate(n=n(),nNoDat=sum(noData)) %>% #All have 205 locations
  filter(nNoDat<(n*0.5)) %>% #Remove days where >50% of data are missing
  select(-noData:-nNoDat) %>% ungroup()

#Overall
p1 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=DO),inherit.aes = FALSE)+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p2 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=Salin),inherit.aes = FALSE)+ #Salin
  scale_colour_distiller(type='div',palette = "YlOrBr")
p3 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=Temp),inherit.aes = FALSE)+ #Temp
  scale_colour_distiller(type='div',palette = "YlOrBr")
p <- ggarrange(p1,p2,p3,ncol=1)
ggsave('./figures/analysis1Figs/wDat_overall_surf.png',p,width=8,height=8)

p1 <- ggmap(basemap)+geom_sf(data=bottomWDat,aes(col=DO),inherit.aes = FALSE)+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p2 <- ggmap(basemap)+geom_sf(data=bottomWDat,aes(col=Salin),inherit.aes = FALSE)+ #Salin
  scale_colour_distiller(type='div',palette = "YlOrBr")
p3 <- ggmap(basemap)+geom_sf(data=bottomWDat,aes(col=Temp),inherit.aes = FALSE)+ #Temp
  scale_colour_distiller(type='div',palette = "YlOrBr")
p <- ggarrange(p1,p2,p3,ncol=1)
ggsave('./figures/analysis1Figs/wDat_overall_bottom.png',p,width=8,height=8)

#Split up maps by date
p <- ggplot()+
  geom_sf(data=surfWDat, aes(col=DO),inherit.aes = FALSE)+
  facet_wrap(~Date)+
  scale_colour_distiller(type='seq',palette = "YlOrBr")
ggsave('./figures/analysis1Figs/wDat_DO_date.png',p,width=8,height=6)

#Time series plots - may have to restrict analysis to specific time "chunks"
p <- surfWDat %>% st_drop_geometry() %>% pivot_longer(c(DO:Salin)) %>% 
  ggplot(aes(x=Date,y=value))+geom_point()+
  facet_wrap(~name,ncol=1)
ggsave('./figures/analysis1Figs/wDat_overall_ts.png',p,width=8,height=8)

png(file = './figures/analysis1Figs/wDat_overall_cor.png',width=6,height=6,units='in',res=150)
surfWDat %>% st_drop_geometry() %>% select(DO:Salin) %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

#How do water values change with depth?
# Salin increases with depth (mixing at surface near rivers), temp decreases (vertical mixing)
# DO decreases with depth, but highly variable in general
p <- wDat %>% st_drop_geometry() %>% 
  pivot_longer(c(DO:Salin)) %>% mutate(fDate=dateCut(Date,6)) %>% 
  ggplot()+  geom_point(aes(x=value,y=Depth),alpha=0.3)+ #Uses absolute depth
  facet_grid(fDate~name,scales='free_x')+
  scale_y_reverse()+labs(y='Depth',x='Measurement')
ggsave('./figures/analysis1Figs/wDat_depth_overall.png',p,width=8,height=10)

#Same thing, but with proportion depth
wDat %>% st_drop_geometry() %>% 
  pivot_longer(c(DO:Salin)) %>% mutate(fDate=dateCut(Date,6)) %>% 
  ggplot()+  geom_point(aes(x=value,y=propDepth),alpha=0.3)+ #Uses proportional depth
  facet_grid(fDate~name,scales='free_x')+
  scale_y_reverse()+labs(y='Proportion Depth',x='Measurement')

# Take a look at spectral data ----------------------------
sDat
head(sDat)

p <- sDat %>% select(YEID:poc,PC1:PC4) %>% 
  pivot_longer(c(chlor_a:poc,PC1:PC4)) %>% 
  filter(!is.na(value)) %>% 
  mutate(name=factor(name,levels=c('chlor_a','nflh','poc','PC1','PC2','PC3','PC4'))) %>% 
  ggplot()+geom_point(aes(x=date_img,y=value))+facet_wrap(~name,ncol=2,scales='free_y')+
  geom_vline(xintercept=range(wDat$Date),col='red',linetype='dashed') #Range of water data
ggsave('./figures/analysis1Figs/sDat_overall_ts.png',p,width=8,height=8)

#Make animation of spectral data across time range
gen_anim <- function() { #Generate set of frames
  sDat3 <- sDat2 %>% mutate(across(chlor_a:poc,~as.vector(scale(.x)))) #%>% pivot_longer(chlor_a:poc)
  mapFun <- function(tau) { #Function to make single frame
    s <- sDat3 %>%  filter(date_img == tau) #subset data to time step tau
      # %>% filter(!is.na(value))
    f <- function(dat,channel,lims){ #Function to make standard figure
      ggmap(basemap) + geom_sf(data=dat,aes(geometry=geometry,colour = {{channel}},size={{channel}}),inherit.aes = FALSE) +
        guides(size=FALSE)+
        scale_size(limits=lims)+
        # scale_colour_gradient(limits=lims)
        scale_colour_distiller(limits=lims,type='div',palette = "YlOrBr",na.value=NA)
    }
    p1 <- f(s,chlor_a,range(sDat3$chlor_a,na.rm=TRUE))
    p2 <- f(s,nflh,range(sDat3$nflh,na.rm=TRUE))
    p3 <- f(s,poc,range(sDat3$poc,na.rm=TRUE))
    p <- ggarrange(p1,p2,p3,ncol=1)
    p <- annotate_figure(p,top = text_grob(as.character(tau),size=10))
    return(p)
  }
  mapFun(tau=sDat2$date_img[1])

  for(t in 1:length(unique(sDat2$date_img))){  # for each time point
    print(mapFun(tau=unique(sDat2$date_img)[t]))   # plot data at this time point
  }
}

setwd("~/Documents/hypoxiaMapping/figures/analysis1Figs")
saveGIF(gen_anim(),movie.name='spectral_anim.gif',interval = 0.5,ani.width=600,ani.height=900)
setwd("~/Documents/hypoxiaMapping")

#Same, but for PCA axes
gen_anim <- function() { #Generate set of frames
  mapFun <- function(tau) { #Function to make single frame
    s <- sDat %>% filter(date_img == tau) #subset data to time step tau
    f <- function(dat,channel,lims){ #Function to make standard figure
      ggmap(basemap) + geom_sf(data=dat,aes(geometry=geometry,colour = {{channel}},size={{channel}}),inherit.aes = FALSE) +
        guides(size=FALSE)+
        scale_size(limits=lims)+
        scale_colour_distiller(limits=lims,type='div',palette = "YlOrBr",na.value=NA)
    }
    p1 <- f(s,PC1,range(sDat$PC1,na.rm=TRUE))
    p2 <- f(s,PC2,range(sDat$PC2,na.rm=TRUE))
    p3 <- f(s,PC3,range(sDat$PC3,na.rm=TRUE))
    p4 <- f(s,PC4,range(sDat$PC4,na.rm=TRUE))
    p <- ggarrange(p1,p2,p3,p4,ncol=1)
    p <- annotate_figure(p,top = text_grob(as.character(tau),size=10))
    return(p)
  }
  mapFun(tau=sDat2$date_img[1])
  
  for(t in 1:length(unique(sDat2$date_img))){  # for each time point
    print(mapFun(tau=unique(sDat2$date_img)[t]))   # plot data at this time point
  }
}
setwd("~/Documents/hypoxiaMapping/figures/analysis1Figs")
saveGIF(gen_anim(),movie.name='spectral_anim_PCA.gif',interval = 0.5,ani.width=600,ani.height=900)
setwd("~/Documents/hypoxiaMapping")

#Hard to see distinct patterns using animations, but it looks like there are "pulses" here and there

#Hovmoller plots for each channel
p1 <- sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% 
  mutate(lon=midcut(lon,20)) %>% st_drop_geometry() %>% 
  mutate(across(chlor_a:poc,~scale(.x))) %>% 
  group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% 
  mutate(lon=as.numeric(as.character(lon))) %>% 
  pivot_longer(chlor_a:poc) %>% 
  ggplot(aes(x=lon,y=date_img))+geom_tile(aes(fill=value))+facet_wrap(~name)+
  scale_fill_distiller(type='div',palette = "YlOrBr",na.value=NA)+
  theme(panel.background = element_rect(fill = 'grey'))+
  labs(x='Longitude',y='Date')+theme(legend.position = 'bottom')

p2 <- sDat2 %>% mutate(lat=st_coordinates(.)[,2]) %>% 
  mutate(lat=midcut(lat,20)) %>% st_drop_geometry() %>% 
  mutate(across(chlor_a:poc,~scale(.x))) %>% 
  group_by(date_img,lat) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% 
  mutate(lat=as.numeric(as.character(lat))) %>% 
  pivot_longer(chlor_a:poc) %>% 
  ggplot(aes(x=lat,y=date_img))+geom_tile(aes(fill=value))+facet_wrap(~name)+
  scale_fill_distiller(type='div',palette = "YlOrBr",na.value=NA)+
  theme(panel.background = element_rect(fill = 'grey'))+
  labs(x='Latitude',y='Date')+theme(legend.position = 'bottom')
(p <- ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend='right'))
ggsave('./figures/analysis1Figs/sDat_hovmoller.png',p,width=12,height=8)

# sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% mutate(lon=cut(lon,breaks=20)) %>% st_drop_geometry() %>% 
#   group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% pull(lon) %>% 
#   levels(.) %>% str()

#Pretty low correlation between classes. Do these actually mean the same thing? Ask YL or LN how these are typically used.
png(file = './figures/analysis1Figs/sDat_overall_cor.png',width=6,height=6,units='in',res=150)
sDat %>% select(chlor_a:poc) %>% na.omit() %>% st_drop_geometry() %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

#Almost no chlor_a measurements on days when DO measurments were taken. Need some way of filling in the gaps
bottomWDat %>% left_join(mutate(st_drop_geometry(sDat),doy=doy),by=c('YEID','doy')) %>% 
  filter(!is.na(chlor_a)) %>% 
  ggplot(aes(chlor_a,DO))+geom_point()

# Fit GAMs of PC1:PC6 ----------------------------------------------------------

#Get smoothers
pcDat <- sDat %>% select(YEID,date_img,E:geometry) %>% filter(!is.na(PC1)) %>% #Strips out missing data
  mutate(sE=sE+rnorm(n(),0,0.1),sN=sN+rnorm(n(),0,0.1)) #Add a tiny bit of noise to make sure that distances between points isn't exactly 0

# "Standard" thin-plate spline approach
# library(parallel)
# cl <- makeCluster(15)
# a <- Sys.time()
# PCmod1 <-bam(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(75,10),d=c(2,1)),data=pcDat,cluster=cl)
# PCmod2 <-bam(PC2~te(sN,sE,doy,bs=c('tp','tp'),k=c(75,10),d=c(2,1)),data=pcDat,cluster=cl)
# PCmod3 <-bam(PC3~te(sN,sE,doy,bs=c('tp','tp'),k=c(75,10),d=c(2,1)),data=pcDat,cluster=cl)
# PCmod4 <-bam(PC4~te(sN,sE,doy,bs=c('tp','tp'),k=c(75,10),d=c(2,1)),data=pcDat,cluster=cl)
# PCmod5 <-bam(PC5~te(sN,sE,doy,bs=c('tp','tp'),k=c(75,10),d=c(2,1)),data=pcDat,cluster=cl)
# PCmod6 <-bam(PC6~te(sN,sE,doy,bs=c('tp','tp'),k=c(75,10),d=c(2,1)),data=pcDat,cluster=cl)
# Sys.time()-a #4.2 mins
# save(pcDat,PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6,file='./data/PCmods.RData')
# par(mfrow=c(2,2)); 
# gam.check(PCmod1); abline(0,1,col='red'); #These look mostly OK
# summary(PCmod1)
# gam.check(PCmod2); abline(0,1,col='red'); #Ehhhh...
# summary(PCmod2)
# gam.check(PCmod3); abline(0,1,col='red'); 
# summary(PCmod3)
# gam.check(PCmod4); abline(0,1,col='red'); #Not great. Looks more like a t-dist. Poor R2
# summary(PCmod4)
# gam.check(PCmod5); abline(0,1,col='red'); 
# summary(PCmod5)
# gam.check(PCmod6); abline(0,1,col='red'); #t-dist again
# summary(PCmod6)
# par(mfrow=c(1,1))

#Approach using soap-film smoothers
shpfileFolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/shapefiles"
regionRect <- st_read(paste0(shpfileFolder,"/region_rectangle.shp"))
coastBuff <- st_read(paste0(shpfileFolder,"/region_coast250kmBuffer.shp")) %>% st_geometry() #250 km buffer zone away from coast
coast <- st_read(paste0(shpfileFolder,"/ne_50m_admin_0_countries_USA.shp")) %>% #Actual coast
  st_crop(regionRect) %>% #Crop to region
  st_geometry()  #Drop other info
dataBoundary <- st_difference(coastBuff,coast)  #Get boundary of analysis area
dataBoundary <- st_sfc(st_polygon(dataBoundary[[1]][[2]][1]),crs=st_crs(dataBoundary)) #Get rid of small islands

# Knots along lines of buffer distances
knotLocs <- mapply(function(d,N){
  dbTemp <- dataBoundary %>% st_transform(3401) %>% st_buffer(-30*1000) %>% st_transform(st_crs(dataBoundary))
  coast[[1]][[1]] %>% st_polygon() %>% st_sfc(.,crs=st_crs(coast)) %>% #Only largest segment
    st_transform(3401) %>% st_buffer(d*1000) %>% st_sf() %>% st_cast('POLYGON') %>%
    mutate(area=as.numeric(st_area(.))) %>% filter(area==max(area)) %>%
    st_geometry() %>% st_transform(st_crs(dataBoundary)) %>%
    st_cast('LINESTRING') %>% st_intersection(dbTemp,.) %>% #Transforms to line, gets intersection with dbTemp
    st_cast('MULTILINESTRING') %>% st_line_merge() %>%  st_transform(3401) %>% #Converts to single line if needed
    st_line_sample(n = N,type='regular') %>% #Samples N points along line
    st_cast('POINT') %>%  st_transform(st_crs(dataBoundary)) #Transforms back to starting crs
},
d=c(50,150,225), #Distances
N=c(25,15,8) #Points within each distance
) %>%
  do.call('c',.)

#Looks OK
plot(dataBoundary)
plot(knotLocs,add=TRUE,pch=19,cex=0.5)
length(knotLocs) #48 knots

#Get knot and boundary locations on same scale
bound <- dataBoundary %>% st_transform(3401) %>%
  st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% dplyr::select(sN,sE)
bound <- list(list(sE=bound$sE,sN=bound$sN)) #Convert to list

kts <- knotLocs %>% st_transform(3401) %>% #Knot locations
  st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% dplyr::select(sE,sN)

library(parallel)
cl <- makeCluster(15)

#Per model: ~10 mins for 60 knots (30 x 2 dims), 30 boundary loops, 15 layers for tensor product
modForm <- "999~te(sN,sE,doy,bs=c('sf','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,30),d=c(2,1))+
  te(sN,sE,doy,bs=c('sw','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,30),d=c(2,1))"
#First te("sf") is for boundary, second te("sw") is for actual soap film

{a <- Sys.time()
PCmod1 <- bam(formula(gsub('999','PC1',modForm)),knots=kts,data=pcDat,cluster=cl) #About 4.3 mins per model
Sys.time()-a}
PCmod2 <- bam(formula(gsub('999','PC2',modForm)),knots=kts,data=pcDat,cluster=cl)
PCmod3 <- bam(formula(gsub('999','PC3',modForm)),knots=kts,data=pcDat,cluster=cl)
PCmod4 <- bam(formula(gsub('999','PC4',modForm)),knots=kts,data=pcDat,cluster=cl)
PCmod5 <- bam(formula(gsub('999','PC5',modForm)),knots=kts,data=pcDat,cluster=cl)
PCmod6 <- bam(formula(gsub('999','PC6',modForm)),knots=kts,data=pcDat,cluster=cl)
stopCluster(cl)
save(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6,file='./data/PCmodsSoap.RData')

par(mfrow=c(2,2));
gam.check(PCmod1); abline(0,1,col='red'); #These look mostly OK
pcDat %>% mutate(resid=resid(PCmod1))

summary(PCmod1)
gam.check(PCmod2); abline(0,1,col='red'); #Ehhhh...
summary(PCmod2)
gam.check(PCmod3); abline(0,1,col='red');
summary(PCmod3)
gam.check(PCmod4); abline(0,1,col='red'); #Not great. Looks more like a t-dist. Poor R2
summary(PCmod4)
gam.check(PCmod5); abline(0,1,col='red');
summary(PCmod5)
gam.check(PCmod6); abline(0,1,col='red'); #t-dist again
summary(PCmod6)
par(mfrow=c(1,1))


# load('./data/PCmods.RData') #Load thin plate splines
load('./data/PCmodsSoap.RData') #Load soap film smoother

# #See how predictions look on spatial grid 
# a <- Sys.time()
# predPC <- with(pcDat,expand.grid(doy=seq(min(doy),max(doy),length.out=9),loc=unique(paste(sE,sN,sep='_')))) %>% 
#   separate(loc,c('sE','sN'),sep='_',convert=TRUE) %>% 
#   mutate(predPC1=predict(PCmod1,newdata=.),sePC1=predict(PCmod1,newdata=.,se.fit=TRUE)$se.fit) %>%
#   mutate(predPC2=predict(PCmod2,newdata=.),sePC2=predict(PCmod2,newdata=.,se.fit=TRUE)$se.fit) %>% 
#   mutate(predPC3=predict(PCmod3,newdata=.),sePC3=predict(PCmod3,newdata=.,se.fit=TRUE)$se.fit) %>% 
#   mutate(predPC4=predict(PCmod4,newdata=.),sePC4=predict(PCmod4,newdata=.,se.fit=TRUE)$se.fit) %>%
#   mutate(predPC5=predict(PCmod5,newdata=.),sePC5=predict(PCmod5,newdata=.,se.fit=TRUE)$se.fit) %>%
#   mutate(predPC6=predict(PCmod6,newdata=.),sePC6=predict(PCmod6,newdata=.,se.fit=TRUE)$se.fit) %>% 
#   mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))
# Sys.time()-a #8 mins
# 
# #Predicted values
# p1 <- predPC %>% select(date,sE,sN,contains('pred')) %>% 
#   pivot_longer(contains('pred')) %>% mutate(name=gsub('pred','',name)) %>% 
#   ggplot()+geom_point(aes(x=sE,y=sN,col=value))+
#   facet_grid(name~date)+
#   scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
#   labs(x='E',y='N',col='PC Value',title='Predicted value')+
#   theme(legend.position='bottom')
# 
# #SE of prediction
# p2 <- predPC %>% select(date,sE,sN,contains('sePC')) %>% 
#   pivot_longer(contains('sePC')) %>% mutate(name=gsub('sePC','PC',name)) %>% 
#   ggplot()+geom_point(aes(x=sE,y=sN,col=value))+
#   facet_grid(name~date)+
#   scale_colour_distiller(type='div',palette = "Reds",direction=1)+
#   labs(x='E',y='N',col='SE of PC Value',title='Standard error of Prediction')+
#   theme(legend.position='bottom')
# 
# #Residual plots
# p3 <- pcDat %>% mutate(residPC1=resid(PCmod1),residPC2=resid(PCmod2),residPC3=resid(PCmod3),
#                        residPC4=resid(PCmod4),residPC5=resid(PCmod5),residPC6=resid(PCmod6)) %>% 
#   pivot_longer(contains('residPC')) %>% mutate(name=gsub('residPC','PC',name)) %>% 
#   mutate(date=dateCut(date_img,9)) %>% 
#   ggplot()+geom_point(aes(x=sE,y=sN,col=value,alpha=abs(value),size=abs(value)))+
#   facet_grid(name~date)+
#   scale_colour_distiller(type='div',palette = "Spectral",direction=1)+
#   scale_alpha_continuous(range=c(0.01,0.75))+
#   labs(x='E',y='N',col='Residual',title='Residual plot')+
#   theme(legend.position='bottom')
# 
# ggsave('./figures/analysis1Figs/PCAmod_pred.png',p1,width=8,height=6)
# ggsave('./figures/analysis1Figs/PCAmod_se.png',p2,width=8,height=6)
# ggsave('./figures/analysis1Figs/PCAmod_resid.png',p3,width=8,height=6)

# Fit model of DO to predicted PCs ----------------------------------

#Get predictions of PCs at all locations through the entire season
predPC <- with(pcDat,expand.grid(YEID=unique(YEID),doy=min(doy):max(doy))) %>% 
  left_join(select(locIndex,-doy,-Date),by='YEID') %>% #Join in spatial info
  mutate(PC1=predict(PCmod1,newdata=.),sePC1=predict(PCmod1,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(PC2=predict(PCmod2,newdata=.),sePC2=predict(PCmod2,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(PC3=predict(PCmod3,newdata=.),sePC3=predict(PCmod3,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(PC4=predict(PCmod4,newdata=.),sePC4=predict(PCmod4,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(PC5=predict(PCmod5,newdata=.),sePC5=predict(PCmod5,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(PC6=predict(PCmod6,newdata=.),sePC6=predict(PCmod6,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

lags <- 0:30 #Try 0 to 30 day lags
fitLagMods <- function(i,interaction=FALSE){

  #Left-join version
  sDatTemp <- predPC %>% #Copy of spectral data
    select(YEID,doy,contains('PC')) %>%
    mutate(doy=doy+i) #Add to go back, subtract to go forward

  sWatTemp <- surfWDat %>% #Copies of water data
    left_join(sDatTemp,by=c('YEID','doy'))
  bWatTemp <- bottomWDat %>%
    left_join(sDatTemp,by=c('YEID','doy'))
  
  # #Choose-this version
  # sWatTemp <- surfWDat #Copies of water data
  # bWatTemp <- bottomWDat
  # chooseThis <- predPC$YEID==sWatTemp$YEID & predPC$doy==(sWatTemp$doy-i) #Location of lagged chlor predictions
  # #Join PCA values onto copies of water data
  # sWatTemp <- predPC %>% filter(chooseThis) %>%
  #   select(contains('PC')) %>% bind_cols(sWatTemp)
  # bWatTemp <- predPC %>% filter(chooseThis) %>%
  #   select(contains('PC')) %>% bind_cols(bWatTemp)
  
  if(!interaction){
    #Fit simple linear models use PC1:5
    sMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=sWatTemp)
    bMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=bWatTemp)
  } else {
    #Interaction models
    sMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=sWatTemp)
    bMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=bWatTemp)
  }
  if(any(is.na(coef(sMod)))|length(coef(sMod))>=nrow(sWatTemp)-5){
    warning('NA coefs or more coefs than data. Model fit singular.')
    return(list(surface=NA,bottom=NA))
  }
  return(list(surface=sMod,bottom=bMod))
}

modList1 <- lapply(lags,fitLagMods,interaction=FALSE)
modList2 <- lapply(lags,fitLagMods,interaction=TRUE)

#Get plots of MSE and R-squared
p1 <- data.frame(lag=lags,surface1=sapply(modList1,function(i) mae(i$surface)),
                 bottom1=sapply(modList1,function(i) mae(i$bottom)),
                 surface2=sapply(modList2,function(i) mae(i$surface)),
                 bottom2=sapply(modList2,function(i) mae(i$bottom))) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+facet_wrap(~name)+
  labs(x='Time lag',y='Mean Absolute Error',col='Model\nType')

p2 <- data.frame(lag=lags,surface1=sapply(modList1,function(i) rmse(i$surface)),
                 bottom1=sapply(modList1,function(i) rmse(i$bottom)),
                 surface2=sapply(modList2,function(i) rmse(i$surface)),
                 bottom2=sapply(modList2,function(i) rmse(i$bottom))) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+facet_wrap(~name)+
  labs(x='Time lag',y='Root Mean Squared Error',col='Model\nType')

p3 <- data.frame(lag=lags,
                 surface1=sapply(modList1,function(i) summary(i$surface)$r.squared),
                 bottom1=sapply(modList1,function(i) summary(i$bottom)$r.squared),
                 surface2=sapply(modList2,function(i) summary(i$surface)$r.squared),
                 bottom2=sapply(modList2,function(i) summary(i$bottom)$r.squared)) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+facet_wrap(~name)+
  labs(x='Time lag',y='R-squared',col='Model\nType')

(p <- ggarrange(p1,p2,p3,ncol=1,common.legend=TRUE,legend='bottom') )
ggsave('./figures/lagPCAmod_smoothed.png',p,width=8,height=8)

# Functional regression using PCAs ----------------------------------------

NdayLag <- 30 #30 days in past
NdayForward <- 0 #0 days in future
dayLags <- -NdayForward:NdayLag

#Matrices to store PCA predictions for past 0:30 days
predMat <- matrix(NA,nrow=nrow(bottomWDat),ncol=length(dayLags),
                       dimnames=list(bottomWDat$YEID,gsub('-','m',paste0('lag',dayLags))))
pcaMatList <- list(PCA1=predMat,PCA2=predMat,PCA3=predMat,PCA4=predMat,PCA5=predMat,PCA6=predMat)

for(p in 1:6){ #PCA dimensions
  for(i in 1:nrow(bottomWDat)){ #For each bottom water measurement
    getDays <- bottomWDat$doy[i]:bottomWDat$doy[i]-dayLags #Which days are 0-30 days behind the measurement?
    pcaMatList[[p]][i,] <- predPC[predPC$YEID == bottomWDat$YEID[i] & predPC$doy %in% getDays,paste0('PC',p)]
  }
}

#Data for functional regression
fdat2 <- list(DO_bottom=bottomWDat$DO,#DO_surf=surfWDat$DO,
             dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
             pcaMat1=pcaMatList$PCA1,pcaMat2=pcaMatList$PCA2,
             pcaMat3=pcaMatList$PCA3,pcaMat4=pcaMatList$PCA4,
             pcaMat5=pcaMatList$PCA5,pcaMat6=pcaMatList$PCA6,
             doy=bottomWDat$doy,sE=bottomWDat$sE,sN=bottomWDat$sN,
             maxDepth=bottomWDat$maxDepth)

basisType <- 'cr' #Thin-plate regression splines with extra shrinkage. Cubic splines have higher R2 but have very strange shapes
#Fit FDA models 
bWatMod2 <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType)+s(dayMat,by=pcaMat5,bs=basisType)+s(dayMat,by=pcaMat6,bs=basisType), 
               data=fdat2) #Bottom water
summary(bWatMod2) #R-squared of about 0.60
par(mfrow=c(2,2)); gam.check(bWatMod2); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
plot(bWatMod2,scheme=1,pages=1)

#Use smoothPred to get FR plots from each smoother

pvals <- unname(round(summary(bWatMod2)$s.table[,4],3))
pvals <- ifelse(pvals==0,'<0.001',paste0('=',as.character(pvals)))

p1 <- lapply(1:6,function(i){
  d <- expand.grid(dayMat=0:30,p=1) #Dataframe
  names(d)[2] <- paste0('pcaMat',i) #Change name of by variable
  smoothPred(m=bWatMod2,dat=d,whichSmooth=i)}) %>% 
  set_names(paste0('PCA',1:6,' (p',pvals,')')) %>% 
  bind_rows(.id='PC') %>% select(-contains('pcaMat')) %>% 
  ggplot(aes(x=dayMat))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred))+facet_wrap(~PC)+geom_hline(yintercept=0,col='red',linetype='dashed')+
  labs(x='Day (lag)',y='Effect')

p2 <- data.frame(pred=predict(bWatMod2),actual=fdat2$DO_bottom) %>% 
  ggplot()+geom_point(aes(x=pred,y=actual))+
  geom_abline(intercept = 0, slope = 1)+
  labs(x='Predicted Bottom DO',y='Actual Bottom DO')

(p <- ggarrange(p1,p2,ncol=2))
ggsave('./figures/frPCA_smoothed.png',p,width=10,height=5)
ggsave('./figures/frPCA_smoothed2.png',p1,width=10,height=5)

# Compare models ----------------------------------------------------------

#RMSE
min(sapply(modList1,function(i) rmse(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) rmse(i$bottom))) #Lagged linear regression - PCA with interactions
rmse(bWatMod2) #Functional regression - PCA

#Which days do these occur on?
lags[which.min(sapply(modList1,function(i) rmse(i$bottom)))] #22-day lag
lags[which.min(sapply(modList2,function(i) rmse(i$bottom)))] #5-day lag

#MAE
min(sapply(modList1,function(i) mae(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) mae(i$bottom))) #Lagged linear regression - PCA with interactions
mae(bWatMod2) #Functional regression - PCA

#R2
max(sapply(modList1,function(i) getR2(i$bottom))) #Lagged linear regression - PCA
max(sapply(modList2,function(i) getR2(i$bottom)),na.rm=TRUE) #Lagged linear regression - PCA
summary(bWatMod2)$r.sq #Functional regression - PCA

#df
sapply(modList1,function(i) getDF(i$bottom))[which.min(sapply(modList1,function(i) rmse(i$bottom)))]
sapply(modList2,function(i) getDF(i$bottom))[which.min(sapply(modList2,function(i) rmse(i$bottom)))]
bWatMod2$df.residual
