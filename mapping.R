#CREATE PREDICTED BOTTOM DO MAPS BASED ON SURFACE MEASUREMENTS
#WRITTEN BY SR, SUMMER 2021

#Load everything -------------------------------------------------

library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(raster)
library(mgcv)
library(parallel)
source('helperFunctions.R')

load('./data/all2014_2.Rdata')
rm(bottomWDat,locIndex,sDat,surfWDat,wDat)

storage <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata"
newfolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/combined"
shpfileFolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/shapefiles"

coastBuff <- st_read(paste0(shpfileFolder,"/region_coast250kmBuffer.shp")) #250 km buffer zone away from coast
# st_read(paste0(shpfileFolder,"/region_rectangle.shp")) %>%
# st_read(paste0(shpfileFolder,"/TWAP_rivers.shp")) %>%
  # ggplot()+geom_sf()


#Get file paths
files <- data.frame(paths=dir(storage,recursive=TRUE,full.names = TRUE)) %>% 
  mutate(file=dir(storage,recursive=TRUE)) %>%
  filter(!grepl('combined',paths)) %>% 
  mutate(file=gsub('.*\\/','',file)) %>% mutate(doy=as.numeric(gsub('\\.tif','',gsub('[AT]2014','',file)))) %>%
  dplyr::select(-file) %>% 
  mutate(platform=ifelse(grepl('terr',paths),'terra','aqua')) %>% 
  mutate(doy=ifelse(platform=='terra',doy-1,doy)) %>% #Terra starts 1 day earlier
  pivot_wider(names_from='platform',values_from = paths) %>%
  mutate(combined=gsub('\\/media.*\\/[AT]','',aqua)) %>% 
  mutate(combined=gsub('2014','2014_',combined)) %>% 
  mutate(combined=paste0(newfolder,'/',combined)) %>% 
  mutate(exists=file.exists(combined))

# #Average between aqua and terra
# a <- Sys.time()
# for(i in 1:nrow(files)){
#   if(!files$exists[i]){
#     aquaDat <- brick(files$aqua[i])
#     terraDat <- brick(files$terra[i])
# 
#     newDat <- lapply(1:dim(aquaDat)[3],function(j) calc(stack(aquaDat[[j]],terraDat[[j]]),fun=mean,na.rm=TRUE))
#     newDat <- brick(newDat)
#     names(newDat) <- names(aquaDat)
#     writeRaster(newDat,filename=files$combined[i],format='GTiff')
#     rm(newDat)
#     gc()
#   }
#   print(paste0('Finished file ',i,' of ',nrow(files)))
# }
# b <- Sys.time()
# b-a #Takes about 10 mins

sDat <- lapply(files$combined,function(datPath){
  cNames <- names(brick(files$aqua[1])) #Channel names
  sD <- brick(datPath) #Read data
  names(sD) <- cNames #Add names
  sD2 <- rasterToPoints(sD) #Convert to point dataframe
  nMissing <- apply(sD2[,!grepl('(x|y)',colnames(sD2))],1,function(x) sum(is.na(x))) #Proportion of missing values in each row (cell)
  sD2 <- data.frame(doy=strsplit(x = datPath,split = c("(\\_|\\.)"))[[1]][2],sD2[nMissing<5,]) #Keep points with 4 or more, and add date
  return(sD2)
}) #Takes about 10 seconds
names(sDat) <- sapply(sDat,function(x) x$doy[1])

#Preview day 152
# sDat[[1]] %>% #Chlor_a on day 152
#   filter(x>(-86.75),x<(-86.25),y<30.2,y>29.8) %>%
#   ggplot(aes(x=x,y=y,fill=chlor_a))+
#   geom_raster()

#Decided to omit non-complete cells - lots of completely missing cells already, so this doesn't hurt that much
sDat2 <- do.call('rbind',sDat) %>% na.omit() %>% #Combine into single DF and remove NAs
  mutate(doy=as.numeric(doy)) %>% #Convert to numeric
  mutate(nflh=ifelse(nflh>(1-1e-5),(1-1e-5),nflh)) %>% #Set upper limit of nflh just below 1
  #Rescales negative values be above 0.95*minimum positive value
  mutate(across(chlor_a:Rrs_678,~ifelse(.x<lwrLimits[names(lwrLimits)==cur_column()],lwrLimits[names(lwrLimits)==cur_column()]*0.95,.x))) 
summary(sDat2)

#Only observations
obs <- sDat2 %>% dplyr::select(-doy:-y) %>% 
  as.matrix() %>% log() #Log-transform

#Calculate PCs 1-6

PCs <- ((obs-outer(rep(1,nrow(obs)),pca1$center))/outer(rep(1,nrow(obs)),pca1$scale)) %*% pca1$rotation
PCs <- PCs[,1:6]
sDat2 <- cbind(sDat2,PCs) #Combine PCs with sDat2
# sDat2 %>% st_as_sf(coords=c("x","y")) %>% filter(doy==153) %>% ggplot()+geom_sf(aes(col=PC1))

#Add coordinate system
sDat2 <- sDat2 %>% st_as_sf(coords=c('x','y')) %>% st_set_crs(4326) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>% #Louisiana offshore
  mutate(sE=(E-mean(unique(E)))/1000,sN=(N-mean(unique(N)))/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% 
  st_transform(4326) %>% dplyr::select(-chlor_a:-sst) %>% unite(loc,E,N) %>% 
  mutate(loc=as.numeric(factor(loc)))

withinBuff <- sDat2 %>% st_intersects(.,coastBuff) %>% sapply(.,function(x) length(x)>0) #Points in sDat2 that are outside of the buffer
sDat2 <- sDat2 %>% filter(withinBuff) #Filter out points outside of buffer
locLookup <- sDat2 %>% dplyr::select(loc:geometry) %>% unique() #Lookup table for locations
# ggplot(locLookup)+geom_sf()

rm(obs,sDat,withinBuff); gc()

# Get predictions from lagged linear model --------------------

load('./data/lagLinMod.Rdata')

# #Function to get predictions from lagged linear model
# getLLpred <- function(l,d,m,interpolate=TRUE,daylag=4){ #l = location ID, d = day
#   newdat <- data.frame(sN=locLookup$sN[which(locLookup$loc==l)],
#                        sE=locLookup$sE[which(locLookup$loc==l)],
#                        doy=d)
#   d <- d-daylag #Lagged day to get data from
#   use <- which(sDat2$doy==d & sDat2$loc==l) #Find matching row in sDat2
#   if(length(use)==1){ #If data exists at that location/time, get from sDat2
#     newdat <- cbind(newdat,sDat2[use,grepl('PC',colnames(sDat2))])
#     newdat$geometry <- NULL
#   } else { #If not, get from GAM models
#     if(interpolate==TRUE){
#       newdat$PC1 <- predict(PCmod1, newdat)
#       newdat$PC2 <- predict(PCmod2, newdat)
#       newdat$PC3 <- predict(PCmod3, newdat)
#       newdat$PC4 <- predict(PCmod4, newdat)
#       newdat$PC5 <- predict(PCmod5, newdat)
#       newdat$PC6 <- predict(PCmod6, newdat)
#     } else {
#       predDO <- NA
#       return(data.frame(DO=NA,interpolated=FALSE))
#     }
#   }
#   predDO <- predict(m,newdat)
#   return(data.frame(DO=predDO,interpolated=length(use)==0))
# }
# 
# #Works
# t(sapply(c(182,196,213,227),function(x){
#   getLLpred(24672,x,m1,interpolate=FALSE)
# },simplify=TRUE))

# a <- Sys.time() #Get predictions at all locations - takes 8 mins: could parallelize though
# preds <- mapply(getLLpred,
#        l=locLookup$loc,d=rep(182,nrow(locLookup)),
#        MoreArgs=list(m=m1),SIMPLIFY=TRUE)
# Sys.time()-a

# #Possible range of predictions: June 4 - July 30 (doy:156-211) 
# #Days to predict on:
# predDays <- c(160,175,190,205) #June 9, 24, July 9, 24 = 160,175,190,205
# a <- Sys.time()
# preds <- mcmapply(getLLpred, #Takes ~9 mins
#                   l=rep(locLookup$loc,length(predDays)),d=rep(predDays,each=nrow(locLookup)),
#                   MoreArgs=list(m=m1,interpolate=FALSE),SIMPLIFY=FALSE,mc.cores=10) %>% 
#   bind_rows()
# Sys.time()-a
# 
# (p1 <- bind_rows(locLookup,locLookup,locLookup,locLookup) %>%
#   bind_cols(preds) %>% 
#   mutate(DO=cut(DO,c(min(DO,na.rm=TRUE),1,2,3,4,5,max(DO,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE)) %>% 
#   mutate(doy=rep(predDays,each=nrow(locLookup))) %>% 
#   mutate(doy=as.Date(paste0('2014-',doy),format='%Y-%j')) %>% 
#   mutate(doy=format(doy,format='%B %d')) %>% 
#   filter(!interpolated) %>% 
#   ggplot()+geom_sf(aes(col=DO))+
#   facet_wrap(~doy,ncol=1)+
#   scale_colour_brewer(type='seq',palette = 'RdBu'))
# ggsave(p1,filename = './figures/mapLLpred.png',width=10,height=10)

#Other solution that just uses "predict"

#Days to average over
dayLag <- 13
d <- seq(min(sDat2$doy+dayLag),max(sDat2$doy+dayLag),by=7)
dayLab <- format(as.Date(paste0('2014-',d),format='%Y-%j'),format='%B %d')
dayLab <- paste(dayLab[1:length(dayLab)-1],':',dayLab[2:length(dayLab)])

(p1 <- sDat2 %>% mutate(DO=predict(m1,.),doy=doy+dayLag,doy2=doy) %>% 
  mutate(doy=cut(doy,breaks=d,labels=dayLab,include.lowest=TRUE)) %>% 
  filter(!is.na(doy)) %>% 
  st_drop_geometry() %>% 
  group_by(doy,loc) %>% summarize(DO=min(DO,na.rm=TRUE)) %>% ungroup() %>% 
  mutate(DO=cut(DO,breaks=c(min(DO,na.rm=TRUE),1,2,3,4,5,max(DO,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE)) %>%
  left_join(locLookup,by='loc') %>% 
  ggplot()+geom_sf(aes(col=DO,geometry=geometry),alpha=0.6)+
  facet_wrap(~doy,ncol=2)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Minimum 1-week lagged-linear predictions'))

ggsave(p1,filename = './figures/mapLLpred2.png',width=14,height=8)
  
#Get predictions from FR model ------------------------------

##Fit PC models
# library(parallel)
# cl <- makeCluster(15)
# #Takes about 10 minutes using 10 cores
# bFuns <- c('tp','tp'); kNum <- c(60,10); dNum <- c(2,1)
# a <- Sys.time()
# PCmod1 <- bam(PC1~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod2 <- bam(PC2~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod3 <- bam(PC3~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod4 <- bam(PC4~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod5 <- bam(PC5~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod6 <- bam(PC6~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# Sys.time()-a
# stopCluster(cl)
# save(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6,file='./data/PCmods_mapping.RData')
load('./data/PCmods_mapping.RData')

load('./data/funRegMod.Rdata')

# l <- 34405
# d <- c(200,211)
l <- sort(unique(locLookup$loc)) #All unique locations
d <- c(196,211) #Days 196,211
daylag <- 30 #30-day lag

cl <- makeCluster(15)
a <- Sys.time() #~1.7 mins for all locations on a single day
newDF <- expand_grid(doy=unique(do.call('c',lapply(d,function(x) x-(daylag:0)))),loc=l) %>% 
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locLookup),by='loc')
newDF <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){require(mgcv); predict.gam(x,newdata=N)},N=newDF) %>% 
  set_names(paste0('predPC',1:6)) %>% bind_cols() %>% bind_cols(newDF,.) %>% 
  mutate(PC1=ifelse(is.na(PC1),predPC1,PC1),PC2=ifelse(is.na(PC2),predPC2,PC2),PC3=ifelse(is.na(PC3),predPC3,PC3)) %>%
  mutate(PC4=ifelse(is.na(PC4),predPC4,PC4),PC5=ifelse(is.na(PC5),predPC5,PC5),PC6=ifelse(is.na(PC6),predPC6,PC6)) 
  # dplyr::select(-contains('pred'),-contains('s'))
Sys.time()-a
stopCluster(cl)

datList <- with(expand.grid(l=l,d=d),list(loc=l,doy=d))
datList$dayMat <- outer(rep(1,length(datList$doy)),0:daylag)

pcaMats <- lapply(paste0('PC',1:6),function(x){
  
  #Works for a single day, but needs to work for multiple days
  
  x <- 'PC1'
  temp <- newDF %>% dplyr::select(doy,loc,x) 
  
  datList$loc
  newDF$loc
  
  datList$doy %in% temp$doy
  
  # #Old method
  #   newDF %>% dplyr::select(doy,loc,x) %>% 
  #   mutate(doy=abs(doy-max(doy))) %>% 
  #   arrange(doy,loc) %>% 
  #   mutate(doy=paste0('d',doy)) %>% 
  #   pivot_wider(values_from = x, names_from = doy) %>% 
  #   column_to_rownames('loc') %>% as.matrix()
}) %>% set_names(paste0('pcaMat',1:6))

datList <- c(datList,pcaMats)

#Predictions are slightly out of the range of data
locLookup %>% arrange(loc) %>% 
  mutate(predDO=predict(bWatMod,newdata=datList)) %>% 
  mutate(predDO=cut(predDO,breaks=c(min(predDO,na.rm=TRUE),1,2,3,4,5,max(predDO,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE)) %>%
  ggplot()+
  geom_sf(aes(geometry=geometry,col=predDO),size=0.2)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title=paste0('Functional regression prediction on day ',d))+
  guides(colour=guide_legend(override.aes = list(size=2)))

# sDat2 %>% mutate(DO=predict(m1,.),doy=doy+dayLag,doy2=doy) %>% 
#     mutate(doy=cut(doy,breaks=d,labels=dayLab,include.lowest=TRUE)) %>% 
#     filter(!is.na(doy)) %>% 
#     st_drop_geometry() %>% 
#     group_by(doy,loc) %>% summarize(DO=min(DO,na.rm=TRUE)) %>% ungroup() %>% 
#     mutate(DO=cut(DO,breaks=c(min(DO,na.rm=TRUE),1,2,3,4,5,max(DO,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE)) %>%
#     left_join(locLookup,by='loc') %>% 
#     ggplot()+  geom_sf(aes(col=DO,geometry=geometry),alpha=0.6)+
#     facet_wrap(~doy,ncol=2)+
#     scale_colour_brewer(type='seq',palette = 'RdBu')+
#     labs(title='Minimum 1-week lagged-linear predictions')



## Range of raw spectral vars for original model (sDat)

## sDat %>% dplyr::select(chlor_a:sst) %>% st_drop_geometry() %>% as.list() %>% sapply(.,range,na.rm=TRUE)
#         chlor_a    nflh      poc    Rrs_412    Rrs_443      Rrs_469     Rrs_488     Rrs_531     Rrs_547     Rrs_555   Rrs_645    Rrs_667    Rrs_678    sst
# [1,]  0.1026964 0.00001    31.16 9.5000e-10 9.5000e-10 1.710095e-05 0.000452001 0.000770001 0.000822001 0.000566001 9.500e-10 9.5000e-10 9.5000e-10 12.645
# [2,] 98.1391449 0.99999 12953.40 1.1678e-02 1.3172e-02 1.628600e-02 0.018146001 0.019711001 0.021056001 0.020276001 1.681e-02 1.5234e-02 1.4894e-02 34.385

#Range for mapping model raw vars
#          chlor_a         nflh     poc    Rrs_412    Rrs_443    Rrs_469     Rrs_488     Rrs_531     Rrs_547     Rrs_555    Rrs_645    Rrs_667    Rrs_678    sst
# [1,]  0.02460281 7.499941e-06    15.0 9.5000e-10 9.5000e-10 1.0000e-09 0.000150001 0.000512001 0.000802001 0.000336001 9.5000e-10 9.5000e-10 9.5000e-10 23.430
# [2,] 99.07888031 1.207985e+00 12953.4 2.1314e-02 1.6558e-02 1.9234e-02 0.020298000 0.024560001 0.025778001 0.024768000 2.0126e-02 1.8378e-02 1.7658e-02 35.615

# nflh >1 and sst higher for mapping vars. nflh can be rescaled, but unclear why sst is higher in larger dataset

## Original Vars
# #After log-transformation
# #apply(sDatMat,2,range,na.rm=TRUE)
#        chlor_a          nflh      poc    Rrs_412    Rrs_443    Rrs_469   Rrs_488   Rrs_531   Rrs_547   Rrs_555    Rrs_645    Rrs_667    Rrs_678      sst
# [1,] -2.275978 -11.51293     3.439135 -20.774559 -20.774559 -10.976377 -7.701826 -7.169119 -7.103769 -7.476915 -20.774559 -20.774559 -20.774559 2.583998
# [2,]  4.586386 -1.000005e-05 9.469114  -4.450048  -4.329662  -4.117449 -4.009305 -3.926578 -3.860570 -3.898317  -4.085781  -4.184225  -4.206797 3.537620

# Mapping vars
# After log-transformation
# apply(obs,2,range)
#        chlor_a          nflh      poc    Rrs_412    Rrs_443    Rrs_469   Rrs_488   Rrs_531   Rrs_547   Rrs_555    Rrs_645    Rrs_667    Rrs_678      sst
# [1,] -3.704894 -11.51293     2.708050 -20.774559 -20.774559 -13.121864 -8.166533 -7.293416 -7.128401 -7.794485 -20.774559 -20.774559 -20.774559 3.158382
# [2,]  4.595916 -1.000005e-05 9.469114  -3.848391  -4.100886  -4.075248 -3.947961 -3.754165 -3.695542 -3.700871  -3.905743  -3.996601  -4.036566 3.572767

## Range of predictor PCs from original FR model (bWatMod)
#           PCA1      PCA2     PCA3      PCA4      PCA5      PCA6
# [1,] -5.951780 -5.278781 -4.41731 -6.213804 -4.474055 -2.894709
# [2,]  9.661089 10.172015  3.28938  2.703307  6.635784  5.546856

## Range of predictors PCs from updated models
## sDat2 %>% st_drop_geometry() %>% dplyr::select(PC1:PC6) %>% as.list() %>% sapply(.,range)
#             PC1       PC2       PC3       PC4        PC5       PC6
# [1,] -0.8211135 -3.103531 -1.599118 -4.384769 -60.526401 -1.172665
# [2,]  0.8562502 10.201641  2.135778  1.345465   8.086439  5.106161

# PC1 = range too small
# PC5 = range too large


#
