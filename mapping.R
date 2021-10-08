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

daylag <- 0 #Uses data from 0 days previous
d <- c(182,189,196,203,211)
l <- sort(unique(locLookup$loc)) #All unique locations
dayLab <- format(as.Date(paste0('2014-',d),format='%Y-%j'),format='%B %d')

newDF <- expand_grid(doy=d,loc=l) %>% #Days/locations to predict at
  mutate(doy2=doy,doy=doy-daylag) %>% #Actual day = doy2, lagged day to get data from = doy
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locLookup),by='loc')

newDF1 <- newDF %>% filter(!is.na(PC1))  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% dplyr::select(-contains('PC')) #DF with NAs

cl <- makeCluster(15)
newDF2 <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){require(mgcv); predict.gam(x,newdata=N)},N=newDF2) %>% 
  set_names(paste0('PC',1:6)) %>% bind_cols() %>% bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN)
stopCluster(cl)
newDF <- bind_rows(newDF1,newDF2) %>% arrange(loc,doy) %>% mutate(predDO=predict(m1,.)) %>% dplyr::select(-sE,-sN) 
rm(newDF1,newDF2) #Cleanup

llMaps <- newDF %>% 
  mutate(predDO=cut(predDO,breaks=c(min(predDO,na.rm=TRUE),1,2,3,4,5,max(predDO,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE)) %>%
  mutate(doy2=factor(doy,labels = dayLab)) %>% 
  left_join(locLookup,by='loc') %>% 
  ggplot()+
    geom_sf(aes(col=predDO,geometry=geometry),size=0.2)+
  facet_wrap(~doy2,ncol=1)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Lagged-linear',col='DO (mg/L)')

  
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


#Possible prediction doy range | lag = days 182 - 211 (Jul 1 - 30)

# "July 01" "July 08" "July 15" "July 22" "July 30"
d <- c(182,189,196,203,211)
# format(as.Date(paste0('2014-',d),format='%Y-%j'),format='%B %d')
l <- sort(unique(locLookup$loc)) #All unique locations
daylag <- 30 #30-day lag

#Dataframe of values to choose from
newDF <- expand_grid(doy=unique(do.call('c',lapply(d,function(x) x-(daylag:0)))),loc=l) %>% 
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locLookup),by='loc')
newDF1 <- newDF %>% filter(!is.na(PC1))  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% dplyr::select(-contains('PC')) #DF with NAs
cl <- makeCluster(15)
a <- Sys.time() #Takes about 2 mins for all 1.1 mil points
newDF2 <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){require(mgcv); predict.gam(x,newdata=N)},N=newDF2) %>% 
  set_names(paste0('PC',1:6)) %>% bind_cols() %>% bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN)
Sys.time()-a
stopCluster(cl)
newDF <- bind_rows(newDF1,newDF2) %>% arrange(loc,doy)
rm(newDF1,newDF2) #Cleanup

#Make list of matrices
datList <- with(expand.grid(l=l,d=d),list(loc=l,doy=d))
datList$dayMat <- outer(rep(1,length(datList$doy)),0:daylag)

nd <- lapply(1:length(d),function(i){
  chooseThese <- (newDF$doy %in% (d[i]-daylag):d[i]) & (newDF$loc %in% l)
  newDF[chooseThese,] %>% data.frame()
}) %>% bind_rows()

pcaMats <- lapply(which(grepl('PC',names(nd))),function(j){
  return(matrix(nd[,j],ncol=ncol(datList$dayMat),byrow=TRUE)[,ncol(datList$dayMat):1])
}) %>% set_names(paste0('pcaMat',1:6))

datList <- c(datList,pcaMats)  
 
#Predictions are slightly out of the range of data
frMaps <- with(datList,data.frame(loc=loc,doy=doy,predDO=predict(bWatMod,newdata=datList))) %>% 
  mutate(predDO=cut(predDO,breaks=c(min(predDO,na.rm=TRUE),1,2,3,4,5,max(predDO,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE)) %>%
  mutate(day=format(as.Date(paste0('2014-',doy),format='%Y-%j'),format='%B %d')) %>% 
  left_join(locLookup,by='loc') %>% 
  ggplot()+
  geom_sf(aes(geometry=geometry,col=predDO),size=0.2)+
  facet_wrap(~day,ncol=1)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Functional regression',col='DO (mg/L)')+
  guides(colour=guide_legend(override.aes = list(size=2)))

library(ggpubr)
p <- ggarrange(llMaps,frMaps,ncol=2,common.legend = TRUE,legend='bottom')

ggsave(p,filename = './figures/mapBothPred.png',width=14,height=8)
