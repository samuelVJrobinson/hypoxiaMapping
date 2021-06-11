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
library(missMDA)
library(gstat)

setwd("~/Documents/hypoxiaMapping")

source('helperFunctions.R')

# Load data from scratch --------------------
#Water data
wDat <- read.csv('./data/sample_waterData.csv') %>% mutate(Date=as.Date(Date,'%Y-%m-%d')) %>%
  mutate(doy=as.numeric(format(Date,format='%j'))) %>%
  # mutate(YEID=ifelse(YEID=='2014_149','2014_003',YEID)) %>% #Set 2014_149 as 2014_003 (see below)
  group_by(YEID) %>% mutate(maxDepth=max(Depth,Depth_dem)) %>% mutate(propDepth=Depth/maxDepth) %>%
  select(-Depth_dem) %>% relocate(YEID:Date,doy,Depth,maxDepth,propDepth,DO,Temp:lon) %>% ungroup() %>%
  st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>% #Louisiana offshore
  mutate(sE=(E-mean(E))/1000,sN=(N-mean(N))/1000) %>% #Center E/N and convert to km
  st_transform(4326)

# #All YEIDs are from a single day
# wDat %>% st_drop_geometry() %>% group_by(YEID) %>% summarize(nDays=length(unique(Date))) %>% filter(nDays!=1)

# # 2014_003, 2014_149 are different times at the same location
# wDat %>% mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2]) %>% unite(loc,x,y) %>% st_drop_geometry() %>% select(YEID,loc) %>%
#   distinct() %>% group_by(loc) %>% mutate(n=n()) %>% filter(n>1)

#Get locations from wDat to join onto spatial data
locIndex <- wDat %>% select(YEID,Date,doy,E:sN,geometry) %>% distinct()

#Spectral data
sDat <- read.csv('./data/sample_spectralData.csv') %>%
  left_join(locIndex,by='YEID') %>% #Join in spatial info
  select(-Date,-doy) %>%
  st_as_sf() %>%
  mutate(date_img=as.Date(date_img,'%Y-%m-%d'),doy=as.numeric(format(date_img,'%j'))) %>%
  # mutate(across(chlor_a:poc,fixNeg2)) %>%  #Change negative values to 95% of smallest + values
  # mutate(nflh=ifelse(nflh>1,1,nflh)) #Limits nflh to less than 1
  mutate(poc=ifelse(poc>quantile(poc,0.008,na.rm=TRUE),poc,quantile(poc,0.008,na.rm=TRUE))) %>%
  mutate(nflh=rescale(fixNeg(nflh,0.95),1e-5,(1-1e-5))) %>% #Rescales nflh to between 0 and 1
  mutate(across(chlor_a:Rrs_678,~rescale(.x,lwr=min(.x[.x>0],na.rm=TRUE)*0.99,upr=NA))) %>% #Rescales values to above 0.99*min + value
  mutate(numNA=apply(st_drop_geometry(.)[,3:16],1,function(x) sum(is.na(x)))) %>% #Count NAs in data columns
  mutate(propNA=numNA/max(numNA))

#choose only surface water, using order of depth measurements
surfWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth)) %>% ungroup() %>%
  filter(depthOrd==1) %>% select(-depthOrd)

#choose only bottom water
bottomWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth,decreasing=TRUE)) %>% ungroup() %>%
  filter(depthOrd==1) %>% select(-depthOrd)

#PCA imputation of missing data imputation for missing data
# Examples: https://link.springer.com/article/10.1007/s11258-014-0406-z
# http://juliejosse.com/wp-content/uploads/2018/05/DataAnalysisMissingR.html

#Most measurements are missing 13 or 14 measurements 
table(sDat$numNA)

#"Most missing" values is chlor_a. Rrs_ values are all missing together
apply(st_drop_geometry(sDat)[,3:16],2,function(x) sum(is.na(x)))/nrow(sDat)

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
# abline(h=0.95,col='red') #Looks like about 5 dims are needed for 95% of variance
# #Singular value decomposition - identical
# svdMat <- svd(scale(sDatMat_noNA)) #First 2 dimensions
# plot(svdMat$d^2/sum(svdMat$d^2),ylab='Relative variance',xlab='Eigenvalue',pch=19)
# plot(cumsum(svdMat$d^2)/sum(svdMat$d^2),ylab='Cumulative variance',xlab='Eigenvalue',pch=19,type='b')
# abline(h=0.95,col='red') #Looks like about 5 dims are needed for 95% of variance
# plot(log(pca$sdev),log(svdMat$d)) #Similar things

# n_components <- estim_ncpPCA(sDatMat, verbose = TRUE,method="Regularized",method.cv="gcv",ncp.min=1,ncp.max=13) #Takes a minute or so
# plot(1:13,n_components$criterion,xlab='Number of Dimensions',ylab='GCV Criterion',type='b',pch=19) #Looks like about 7 components is OK for prediction
sDatMat_imputed <- imputePCA(sDatMat,ncp=7,scale=TRUE,method='Regularized') #Impute missing data using 7 dimensions 

head(sDatMat_imputed$completeObs) #Filled-in data
head(sDatMat) #Original data

pca1 <- prcomp(sDatMat_imputed$completeObs,scale=TRUE) #Get PCA values
plot(1:length(pca1$sdev),cumsum(pca1$sdev^2)/sum(pca1$sdev^2),type='b',xlab='Principle Component',ylab='Cumulative Variance',pch=19)
abline(h=0.95,col='red') #Looks like about 4 dims are needed for 95% of variance, 5 for 97.5%
abline(h=0.975,col='red')

#Factor loadings of first 5 PCs
pca1$rotation[,1:5] %>% data.frame() %>% rownames_to_column(var='var') %>%
  pivot_longer(PC1:PC5) %>%
  mutate(name=factor(name,labels=paste0('PC',1:5,': ',round(pca1$sdev[1:5]^2/sum(pca1$sdev^2),3)*100,'% Variance'))) %>%
  mutate(var=factor(var,levels=rev(unique(var)))) %>%
  ggplot()+
  # geom_point(aes(y=var,x=value))+
  geom_col(aes(y=var,x=value))+geom_vline(xintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~name)+
  labs(x='Loading',y=NULL,title='Variable loadings for Principle Components 1-5')

#Join imputed PCA values back to original dataset (NA if no data on that day)
sDat_pca <- pca1$x[,1:5] %>% as.data.frame() %>% rownames_to_column('ID')
sDat <- sDat %>% unite(ID,YEID,date_img,sep='_',remove=FALSE) %>% 
  left_join(sDat_pca,by='ID') %>% select(-ID) %>% mutate(gap=is.na(PC1))
rm(sDat_pca,sDatMat_imputed,sDatMat) #cleanup
save(bottomWDat,locIndex,pca1,sDat,surfWDat,wDat,file='./data/all2014.Rdata')

#Load data from saved file
load('./data/all2014.Rdata')

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

# Fit GAM of chlor_a ----------------------------------------------------

#Using entire dataset. Fit is somewhat better if using reduced dataset (sDat2), but going with this for now.
#Would fit be better using reduced dataset because of 
chlorDat <- sDat %>% select(YEID:chlor_a,E:geometry) %>% 
  mutate(doy=format(date_img,format='%j'),doy=as.numeric(doy)) %>% 
  mutate(logChlor=log(chlor_a)) %>% 
  filter(!is.na(chlor_a)) #Strips out missing data

# library(parallel) #No significant difference in fitting times using bam
# detectCores()
# cl <- makeCluster(12)
# stopCluster(cl) 

#Uses 50 bases for N,E and 10 for time. Tensor product = 500 terms
chlorMod <- gam(logChlor~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=chlorDat)
# chlorMod <- gam(chlor_a~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=chlorDat,family=Gamma) #Takes a ton of time
# chlorMod2 <- gam(chlor_a~te(sN,sE,doy,bs=c('tp','tp'),k=c(30,5),d=c(2,1)),data=chlorDat,family=Gamma())
# chlorMod3 <- gam(chlor_a~te(sN,sE,doy,bs=c('tp','tp'),k=c(30,5),d=c(2,1)),data=chlorDat,family=inverse.gaussian()) 
# par(mfrow=c(2,2)); gam.check(chlorMod2); abline(0,1,col='red'); par(mfrow=c(1,1))
# par(mfrow=c(2,2)); gam.check(chlorMod3); abline(0,1,col='red'); par(mfrow=c(1,1))

summary(chlorMod)
par(mfrow=c(2,2)); gam.check(chlorMod); abline(0,1,col='red'); par(mfrow=c(1,1))
plot(chlorMod,scheme=3,too.far=0.1)
chlorDat$resid <- resid(chlorMod) #Residuals

#Prediction grid
predChlor <- with(chlorDat,expand.grid(doy=seq(min(doy),max(doy),length.out=9),loc=unique(paste(sE,sN,sep='_')))) %>% 
  separate(loc,c('sE','sN'),sep='_',convert=TRUE) %>% 
  mutate(pred=predict(chlorMod,newdata=.),se=predict(chlorMod,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

p1 <- predChlor %>% ggplot()+geom_point(aes(x=sE,y=sN,col=pred))+
  facet_wrap(~date,nrow=1)+
  scale_colour_distiller(type='div',palette = "Greens",direction=1)+
  labs(x='E',y='N',col='Predicted\nlog(Chlor)')+
  theme(legend.position='bottom')

p2 <- predChlor %>% ggplot()+geom_point(aes(x=sE,y=sN,col=se))+
  facet_wrap(~date,nrow=1)+
  scale_colour_distiller(type='div',palette = "YlOrRd",direction=1) +
  labs(x='E',y='N',col='SE')+
  theme(legend.position='bottom')

(p <- ggarrange(p1,p2,nrow=2))
ggsave('./figures/analysis1Figs/chlorMod.png',p,width=11,height=6)

#Plot of residuals over space for 12 time slices
chlorDat %>% mutate(cDate=dateCut(date_img,12)) %>%
  ggplot()+geom_sf(aes(geometry=geometry,col=resid,size=abs(resid)),alpha=0.5)+
  facet_wrap(~cDate)+ scale_colour_distiller(type='div',palette = "RdBu")

#Variogram of residuals

#Overall spatial autocorrelation
p1 <- select(chlorDat,YEID,date_img,resid) %>% 
  st_transform(3401) %>% as('Spatial') %>% 
  variogram(resid~1, data=.,width=2000) %>% 
  ggplot()+geom_line(aes(x=dist,y=gamma))+
  labs(x='Distance',y='Semivariance',title='Overall')
#First day
p2 <- chlorDat %>% filter(date_img==first(date_img)) %>% 
  select(YEID,date_img,resid) %>% 
  st_transform(3401) %>% as('Spatial') %>% 
  variogram(resid~1, data=.,width=2000) %>% 
  ggplot()+geom_line(aes(x=dist,y=gamma))+
  labs(x='Distance',y='Semivariance',title='First Day')
#Last day
p3 <- chlorDat %>% filter(date_img==last(date_img)) %>% 
  select(YEID,date_img,resid) %>% 
  st_transform(3401) %>% as('Spatial') %>% 
  variogram(resid~1, data=.,width=2000) %>% 
  ggplot()+geom_line(aes(x=dist,y=gamma))+
  labs(x='Distance',y='Semivariance',title='Last Day')
ggarrange(p1,p2,p3,nrow=3)

#Spatial autocorr (Moran's I) at each day - strong spatial autocorrelation
m <- sapply(sort(unique(chlorDat$date_img)),function(i){
  temp <- chlorDat %>% filter(date_img==i)
  
  if(nrow(temp)<=5) return(NA)
  
  d <- temp %>% st_distance()
  units(d) <- NULL
  d <- 1000*d/max(d) #Change to km
  d <- 1/d #inverse
  diag(d) <- 0
  r <- temp %>% pull(resid) #Residuals
  if(sum(!is.finite(d))>0){ #IF NA/Inf (sites at same location)
    notFinite <- which(apply(d,2,function(x) as.logical(sum(!is.finite(x)))))[-1]
    d <- d[-notFinite,-notFinite] #Remove rows with NA/Inf
    r <- r[-notFinite]
  }
  if(!is.finite(ape::Moran.I(r,d)$p.value)) stop('NaN value')
  
  return(ape::Moran.I(r,d)$p.value)
})
m <- m[!is.na(m)] #Remove NAs
sum(m<(0.05/length(m)))/length(m) #73% of days show evidence of spatial autocorrelation in residuals

#Temporal autocorrelation (Durbin-Watson) at each location - uses order rather than day
n <- sapply(unique(chlorDat$YEID),function(i){
  temp <- chlorDat %>% filter(YEID==i) %>% lmtest::dwtest(resid~1,data=.)
  return(temp$p.value)
})
sum(n<(0.05/length(n)))/length(n) #Not much

# #ACF of residuals. Creates ts() object with NAs in missing days
# temp <- chlorDat %>% filter(YEID==last(YEID)) %>% st_drop_geometry() %>% mutate(doy=1+doy-min(doy))
# r <- rep(NA,max(temp$doy))
# r[temp$doy] <- temp$resid
# a <- acf(r,na.action = na.pass,type='partial',lag.max=30)
# 
# r <- ts(r,start=1,end=length(r))
# ar(r,na.action='na.pass')
# Box.test(r)

#Empirical ST variogram of residuals
a <- Sys.time()
vv <- chlorDat %>% slice_sample(n=1000) %>% 
  st_drop_geometry() %>% 
  stConstruct(x=.,space=c('E','N'),time='date_img') %>% 
  variogramST(resid~1, #Takes about 2 mins if using 1000 data points, width=2000, tlags= 100, cores =10
            data=.,
            width=10000,
            tlags=seq(0,30,by=1),cores = 10,tunit='days')
Sys.time()-a

vv %>% 
  ggplot(aes(x=spacelag/1000,y=timelag,fill=gamma))+geom_raster() +
  labs(x='Distance (km)',y='Time (days)',fill='Semi-\nvariance')

# #Try larger tensor product (10x10x10).  What has 500 more coefficients gotten us?
# chlorMod2 <- gam(logChlor~te(sN,sE,doy,k=10),data=chlorDat)  #Better AIC, but residual checks are about the same.
# summary(chlorMod)
# summary(chlorMod2)
# anova(chlorMod,chlorMod2,test='Chisq') #Slightly worse CV (0.347 vs 0.350)
# AIC(chlorMod,chlorMod2) #AIC is worse
# chlorDat$resid2 <- resid(chlorMod2)
# 
# #Plot of residuals over space for 5 time slices. Looks like both models are doing poorly at the same places
# pivot_longer(chlorDat,resid:resid2) %>%
#   mutate(cDate=dateCut(date_img,6)) %>%
#   ggplot()+geom_point(aes(x=E,y=N,col=value,size=abs(value)),alpha=0.5)+
#   facet_grid(cDate~name)+
#   scale_colour_distiller(type='div',palette = "RdBu")
# 
# #Mod1 is worse around the edges of the distribution, while mod2 is worse at the centre
# chlorDat %>% mutate(resDiff=abs(resid)-abs(resid2)) %>%  #Difference of absolute residuals
#   mutate(cDate=dateCut(date_img,12)) %>%
#   ggplot()+geom_point(aes(x=E,y=N,col=resDiff,size=abs(resDiff)),alpha=0.5)+
#   facet_wrap(~cDate)+
#   scale_colour_distiller(type='div',palette = "RdBu") #Red = model1 worse, Blue = model2 worse
# 
# #Looks about the same in time
# pivot_longer(chlorDat,resid:resid2) %>%
#   ggplot(aes(x=date_img,y=value))+geom_point()+
#   facet_wrap(~name,ncol=1)

# Predict DO using lagged chlor_a ------------------------------------------------

#Simple model using lagged chlor_a to predict DO for both 

#Get predictions at all locations through the entire season
predChlor <- with(chlorDat,expand.grid(YEID=unique(YEID),doy=min(doy):max(doy))) %>% 
  left_join(select(locIndex,-doy,-Date),by='YEID') %>% #Join in spatial info
  mutate(pred=predict(chlorMod,newdata=.),se=predict(chlorMod,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

lags <- 0:30 #Try 0 to 40 day lags (negative = days after, positive = days before)
modList <- lapply(lags,function(i){
  #Copies of water data
  sWatTemp <- surfWDat
  bWatTemp <- bottomWDat 
  chooseThis <- predChlor$YEID==sWatTemp$YEID & predChlor$doy==(sWatTemp$doy-i) #Location of lagged chlor predictions
  sWatTemp$lagChlor <- predChlor$pred[chooseThis]
  sWatTemp$lagDay <- predChlor$doy[chooseThis]
  bWatTemp$lagChlor <- predChlor$pred[chooseThis]
  bWatTemp$lagDay <- predChlor$doy[chooseThis]
  sMod <- lm(DO~lagChlor,data=sWatTemp)
  bMod <- lm(DO~lagChlor,data=bWatTemp)
  return(list(surface=sMod,bottom=bMod))
})

p1 <- data.frame(lag=lags,surface=sapply(lapply(modList,function(x) x$surface),rmse),bottom=sapply(lapply(modList,function(x) x$bottom),rmse)) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='RMSE')

p2 <- data.frame(lag=lags,surface=sapply(lapply(modList,function(x) x$surface),mae),bottom=sapply(lapply(modList,function(x) x$bottom),mae)) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='MAE')

p3 <- data.frame(lag=lags,surface=sapply(modList,function(i) summary(i$surface)$r.squared),
           bottom=sapply(modList,function(i) summary(i$bottom)$r.squared)) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='R-squared')

ggarrange(p1,p2,p3,ncol=1)

#Take a look at model on day -11
par(mfrow=c(2,1))
dayLag <- which.min(sapply(modList,function(i) sum(resid(i$bottom)^2)))
bestMod <- modList[[dayLag]]$bottom #Bottom water model
with(bestMod$model,plot(lagChlor,DO,xlab=paste0('Chlor_a at day ',lags[dayLag]),ylab='Bottom DO'))
abline(bestMod) #Fairly good relationship (not exactly linear, but OK for now)

dayLag <- which.min(sapply(modList,function(i) sum(resid(i$surface)^2)))
bestMod <- modList[[dayLag]]$surface #Surface water model
with(bestMod$model,plot(lagChlor,DO,xlab=paste0('Chlor_a at day ',lags[dayLag]),ylab='Surface DO'))
abline(bestMod) #No real relationship
par(mfrow=c(1,1))

# NOTES:
# Looks like surface water DO is more poorly predicted by chlor_a. However, lower MSE, so probably just less variable overall
# Strangely, future chlor_a is a slightly better predictor than past chlor_a. What is going on here?
# Hypotheses:
# 1) Just a weird dataset, and won't show up in future years
# 2) Some thing is actually causing a DO drop, followed by an increase in chlor_a

#Predict DO using FR of chlor_a -----------------------------------------------

NdayLag <- 30 #30 days in past
NdayForward <- 0 #20 days in future
dayLags <- -NdayForward:NdayLag

#Matrix to store chlor_a predictions for past 0:30 days
predChlorMat <- matrix(NA,nrow=nrow(bottomWDat),ncol=length(dayLags),
                       dimnames=list(bottomWDat$YEID,gsub('-','m',paste0('lag',dayLags))))

for(i in 1:nrow(bottomWDat)){
  getDays <- bottomWDat$doy[i]:bottomWDat$doy[i]-dayLags
  predChlorMat[i,] <- predChlor$pred[predChlor$YEID == bottomWDat$YEID[i] & predChlor$doy %in% getDays]
}

#Data for functional regression
fdat <- list(DO_bottom=bottomWDat$DO,DO_surf=surfWDat$DO,
          dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
          chlorMat=predChlorMat,
          doy=bottomWDat$doy,sE=bottomWDat$sE,sN=bottomWDat$sN,
          maxDepth=bottomWDat$maxDepth)

#Fit FDA models - looks slightly better, but still only gets R2 of 0.3, so not much better than simple lag model
bWatMod <- gam(DO_bottom ~ s(dayMat,by=chlorMat), data=fdat) #Bottom water
summary(bWatMod)
par(mfrow=c(2,2)); gam.check(bWatMod); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
plot(bWatMod,scheme=1,pages=1); abline(h=0,col='red')

data.frame(pred=predict(bWatMod),actual=fdat$DO_bottom) %>% 
  ggplot()+geom_point(aes(x=pred,y=actual))+
  geom_abline(intercept = 0, slope = 1)+
  labs(x='Predicted Bottom DO',y='Actual Bottom DO')

p1 <- bottomWDat %>% mutate(resid=resid(bWatMod)) %>% st_jitter(factor=0.01) %>% 
  ggplot()+geom_sf(aes(geometry=geometry,col=resid,size=abs(resid)),alpha=0.5)+
  scale_colour_distiller(type='div',palette = "RdBu")+theme(legend.position='bottom')+
  guides(size=FALSE)
p2 <- bottomWDat %>% mutate(resid=resid(bWatMod)) %>% ggplot()+geom_point(aes(x=doy,y=resid))+geom_hline(yintercept = 0)
ggarrange(p1,p2,ncol=1)

sWatMod <- gam(DO_surf ~ s(dayMat,by=chlorMat), data=fdat) #Surface water DO
summary(sWatMod)
par(mfrow=c(2,2)); gam.check(sWatMod); abline(0,1,col='red'); par(mfrow=c(1,1))
plot(sWatMod,pages=1,scheme=1); abline(h=0,col='red') 

p1 <- surfWDat %>% mutate(resid=resid(sWatMod)) %>% st_jitter(factor=0.01) %>% 
  ggplot()+geom_sf(aes(geometry=geometry,col=resid,size=abs(resid)),alpha=0.5)+
  scale_colour_distiller(type='div',palette = "RdBu")+theme(legend.position='bottom')+
  guides(size=FALSE)
p2 <- surfWDat %>% mutate(resid=resid(sWatMod)) %>% ggplot()+geom_point(aes(x=doy,y=resid))+geom_hline(yintercept = 0)
ggarrange(p1,p2,ncol=1)

# Fit GAMs of PC1:PC4 ----------------------------------------------------------

# load('./data/PCmods.RData')

pcDat <- sDat %>% select(YEID,date_img,E:geometry) %>% filter(!is.na(PC1)) %>% #Strips out missing data
  mutate(sE=sE+rnorm(n(),0,0.1),sN=sN+rnorm(n(),0,0.1)) #Add a bit of random noise to make sure that distances between points aren't 0

PCmod1 <-gam(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat)
PCmod2 <-gam(PC2~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat)
PCmod3 <-gam(PC3~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat) 
PCmod4 <-gam(PC4~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat) 
PCmod5 <-gam(PC5~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat) 
save(pcDat,PCmod1,PCmod2,PCmod3,PCmod4,file='./data/PCmods.RData')

par(mfrow=c(2,2)); 
gam.check(PCmod1); abline(0,1,col='red'); #These look mostly OK
gam.check(PCmod2); abline(0,1,col='red'); 
gam.check(PCmod3); abline(0,1,col='red'); 
gam.check(PCmod4); abline(0,1,col='red'); #Bunch of lower outliers
gam.check(PCmod5); abline(0,1,col='red'); 
par(mfrow=c(1,1))

# #Haven't gotten this to complete fitting yet. Also doesn't account for spatial correlation
# PCmod1 <-gamm(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat,
#               correlation=corExp(form=~sN+sE,nugget=TRUE))

summary(PCmod1) 
summary(PCmod2) 
summary(PCmod3) 
summary(PCmod4) #Not as good
summary(PCmod5) #Not as good

#See how predictions look on spatial grid - takes a few minutes
predPC <- with(pcDat,expand.grid(doy=seq(min(doy),max(doy),length.out=9),loc=unique(paste(sE,sN,sep='_')))) %>% 
  separate(loc,c('sE','sN'),sep='_',convert=TRUE) %>% 
  mutate(predPC1=predict(PCmod1,newdata=.),sePC1=predict(PCmod1,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(predPC2=predict(PCmod2,newdata=.),sePC2=predict(PCmod2,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(predPC3=predict(PCmod3,newdata=.),sePC3=predict(PCmod3,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(predPC4=predict(PCmod4,newdata=.),sePC4=predict(PCmod4,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(predPC5=predict(PCmod5,newdata=.),sePC5=predict(PCmod5,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

#Predicted values
p1 <- predPC %>% select(date,sE,sN,contains('pred')) %>% 
  pivot_longer(contains('pred')) %>% mutate(name=gsub('pred','',name)) %>% 
  ggplot()+geom_point(aes(x=sE,y=sN,col=value))+
  facet_grid(name~date)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
  labs(x='E',y='N',col='PC Value',title='Predicted value')+
  theme(legend.position='bottom')

#SE of prediction
p2 <- predPC %>% select(date,sE,sN,contains('sePC')) %>% 
  pivot_longer(contains('sePC')) %>% mutate(name=gsub('sePC','PC',name)) %>% 
  ggplot()+geom_point(aes(x=sE,y=sN,col=value))+
  facet_grid(name~date)+
  scale_colour_distiller(type='div',palette = "Reds",direction=1)+
  labs(x='E',y='N',col='SE of PC Value',title='Standard error of Prediction')+
  theme(legend.position='bottom')

#Residual plots
p3 <- pcDat %>% mutate(residPC1=resid(PCmod1),residPC2=resid(PCmod2),residPC3=resid(PCmod3),residPC4=resid(PCmod4),residPC5=resid(PCmod5)) %>% 
  pivot_longer(contains('residPC')) %>% mutate(name=gsub('residPC','PC',name)) %>% 
  mutate(date=dateCut(date_img,9)) %>% 
  ggplot()+geom_point(aes(x=sE,y=sN,col=value,alpha=abs(value),size=abs(value)))+
  facet_grid(name~date)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=1)+
  scale_alpha_continuous(range=c(0.01,0.75))+
  labs(x='E',y='N',col='Residual',title='Residual plot')+
  theme(legend.position='bottom')

ggsave('./figures/analysis1Figs/PCAmod_pred.png',p1,width=8,height=6)
ggsave('./figures/analysis1Figs/PCAmod_se.png',p2,width=8,height=6)
ggsave('./figures/analysis1Figs/PCAmod_resid.png',p3,width=8,height=6)

#Fit model of DO to predicted PCs ----------------------------------

#Get predictions of PCs at all locations through the entire season
predPC <- with(pcDat,expand.grid(YEID=unique(YEID),doy=min(doy):max(doy))) %>% 
  left_join(select(locIndex,-doy,-Date),by='YEID') %>% #Join in spatial info
  mutate(predPC1=predict(PCmod1,newdata=.),sePC1=predict(PCmod1,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(predPC2=predict(PCmod2,newdata=.),sePC2=predict(PCmod2,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(predPC3=predict(PCmod3,newdata=.),sePC3=predict(PCmod3,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(predPC4=predict(PCmod4,newdata=.),sePC4=predict(PCmod4,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

lags <- 0:30 #Try 0 to 30 day lags (negative = days after, positive = days before)
modList2 <- lapply(lags,function(i){
  sWatTemp <- surfWDat #Copies of water data
  bWatTemp <- bottomWDat 
  chooseThis <- predPC$YEID==sWatTemp$YEID & predPC$doy==(sWatTemp$doy-i) #Location of lagged chlor predictions
  #Join PCA values onto copies of water data
  sWatTemp <- predPC %>% filter(chooseThis) %>% 
    select(contains('predPC')) %>% bind_cols(sWatTemp)
  bWatTemp <- predPC %>% filter(chooseThis) %>% 
    select(contains('predPC')) %>% bind_cols(bWatTemp)
  #Fit simple linear models use PC1:4
  sMod <- lm(DO~predPC1+predPC2+predPC3+predPC4,data=sWatTemp)
  bMod <- lm(DO~predPC1+predPC2+predPC3+predPC4,data=bWatTemp)
  
  # # Slightly better prediction using interactions (R2: 0.47 vs 0.36, MSE: 353 vs 455)
  # sMod <- lm(DO~predPC1*predPC2*predPC3*predPC4,data=sWatTemp)
  # bMod <- lm(DO~predPC1*predPC2*predPC3*predPC4,data=bWatTemp)
  return(list(surface=sMod,bottom=bMod))
})

#Get plots of MSE and R-squared
p1 <- data.frame(lag=lags,surface=sapply(modList2,function(i) mean(abs(resid(i$surface)))),
                 bottom=sapply(modList2,function(i) mean(abs(resid(i$bottom))))) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='Mean Absolute Error')

p2 <- data.frame(lag=lags,surface=sapply(modList2,function(i) summary(i$surface)$r.squared),
                 bottom=sapply(modList2,function(i) summary(i$bottom)$r.squared)) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='R-squared')
ggarrange(p1,p2,ncol=1) 

#Best MSE/R2 for bottom DO 8 days in the past, 25 days in past for surface DO
lags[which.min(sapply(modList2,function(i) sum(resid(i$bottom)^2)))]
lags[which.min(sapply(modList2,function(i) sum(resid(i$surface)^2)))]

#Take a look at both models
par(mfrow=c(2,1))
dayLag <- which.min(sapply(modList2,function(i) sum(resid(i$bottom)^2)))
bestMod <- modList2[[dayLag]]$bottom #Bottom water model
plot(bestMod$model$DO,predict(bestMod),xlab='Actual DO',ylab='Predicted DO',main=paste0('Bottom DO (',lags[dayLag],' day lag)'))
abline(0,1) #Fairly good relationship (not linear, but OK for now)

dayLag <- which.min(sapply(modList2,function(i) sum(resid(i$surface)^2)))
bestMod <- modList2[[dayLag]]$surface #Surface water model
plot(bestMod$model$DO,predict(bestMod),xlab='Actual DO',ylab='Predicted DO',main=paste0('Surface DO (',lags[dayLag],' day lag)'))
abline(0,1) #No real relationship
par(mfrow=c(1,1))

#Get best models
bottomMod <- modList2[[which.min(sapply(modList2,function(i) sum(resid(i$bottom)^2)))]]$bottom #Bottom water model
surfaceMod <- modList2[[which.min(sapply(modList2,function(i) sum(resid(i$surface)^2)))]]$surface #surface water model

summary(bottomMod)
par(mfrow=c(2,1)); plot(bottomMod,which=c(1,2)); par(mfrow=c(1,1)) #Not too bad here
par(mfrow=c(2,1)); plot(surfaceMod,which=c(1,2)); par(mfrow=c(1,1)) #Problems here

# Functional regression using PCAs ----------------------------------------

NdayLag <- 30 #30 days in past
NdayForward <- 0 #0 days in future
dayLags <- -NdayForward:NdayLag

#Matrices to store PCA predictions for past 0:30 days
predMat <- matrix(NA,nrow=nrow(bottomWDat),ncol=length(dayLags),
                       dimnames=list(bottomWDat$YEID,gsub('-','m',paste0('lag',dayLags))))
pcaMatList <- list(PCA1=predMat,PCA2=predMat,PCA3=predMat,PCA4=predMat)

for(p in 1:4){ #PCA dimensions
  for(i in 1:nrow(bottomWDat)){ #For each bottom water measurement
    getDays <- bottomWDat$doy[i]:bottomWDat$doy[i]-dayLags #Which days are 0-30 days behind the measurement?
    pcaMatList[[p]][i,] <- predPC[predPC$YEID == bottomWDat$YEID[i] & predPC$doy %in% getDays,paste0('predPC',p)]
  }
}

#Data for functional regression
fdat2 <- list(DO_bottom=bottomWDat$DO,DO_surf=surfWDat$DO,
             dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
             pcaMat1=pcaMatList$PCA1,pcaMat2=pcaMatList$PCA2,
             pcaMat3=pcaMatList$PCA3,pcaMat4=pcaMatList$PCA4,
             doy=bottomWDat$doy,sE=bottomWDat$sE,sN=bottomWDat$sN,
             maxDepth=bottomWDat$maxDepth)

basisType <- 'ts' #Thin-plate regression splines with extra shrinkage. Cubic splines have higher R2 but have very strange shapes
#Fit FDA models 
bWatMod2 <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType), 
               data=fdat2) #Bottom water
summary(bWatMod2) #R-squared of about 0.46
par(mfrow=c(2,2)); gam.check(bWatMod2); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
plot(bWatMod2,scheme=1,pages=1)

#Use smoothPred to get FR plots from each smoother
p1 <- lapply(1:4,function(i){
  d <- expand.grid(dayMat=0:30,p=1) #Dataframe
  names(d)[2] <- paste0('pcaMat',i) #Change name of by variable
  smoothPred(m=bWatMod2,dat=d,whichSmooth=i) 
}) %>% set_names(paste0('PCA',1:4)) %>% bind_rows(.id='PC') %>% 
  select(-contains('pcaMat')) %>% 
  ggplot(aes(x=dayMat))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred))+facet_wrap(~PC)+geom_hline(yintercept=0,col='red',linetype='dashed')+
  labs(x='Day (lag)',y='Effect')

p2 <- data.frame(pred=predict(bWatMod2),actual=fdat2$DO_bottom) %>% 
  ggplot()+geom_point(aes(x=pred,y=actual))+
  geom_abline(intercept = 0, slope = 1)+
  labs(x='Predicted Bottom DO',y='Actual Bottom DO')

(p <- ggarrange(p1,p2,ncol=2))
ggsave('./figures/analysis1Figs/PCA_FRmod.png',p,width=10,height=5)

#TO DO: 
# See how shape of PC1:PC4 curves changes with time (look up how this is done in arthropod project)


# Compare models ----------------------------------------------------------

min(sapply(modList,function(i) maa(i$bottom))) #Lagged linear regression - chlor_a
maa(bWatMod) #Functional regression - chlor_a
min(sapply(modList2,function(i) maa(i$bottom))) #Lagged linear regression - PCA
maa(bWatMod2) #Functional regression - PCA

#R2
max(sapply(modList,function(i) summary(i$bottom)$adj.r.squared)) #Lagged linear regression - chlor_a
summary(bWatMod)$r.sq #Functional regression - chlor_a
max(sapply(modList2,function(i) summary(i$bottom)$adj.r.squared)) #Lagged linear regression - PCA
summary(bWatMod2)$r.sq #Functional regression - PCA

