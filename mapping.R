#CREATE PREDICTED BOTTOM DO MAPS BASED ON SURFACE MEASUREMENTS
#WRITTEN BY SR, SUMMER 2021

#Load everything -------------------------------------------------

library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(raster)
library(mgcv)
library(parallel)
library(beepr)
library(ggpubr)
source('helperFunctions.R')

load('./data/all2014_2.Rdata')
rm(bottomWDat,locIndex,sDat,surfWDat,wDat)

storage <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata"
newfolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/combined"
shpfileFolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/shapefiles"

regionRect <- st_read(paste0(shpfileFolder,"/region_rectangle.shp"))

coastBuff <- st_read(paste0(shpfileFolder,"/region_coast250kmBuffer.shp")) %>% st_geometry() #250 km buffer zone away from coast
coast <- st_read(paste0(shpfileFolder,"/ne_50m_admin_0_countries_USA.shp")) %>% #Actual coast
  st_crop(regionRect) %>% #Crop to region
  st_geometry()  #Drop other info
dataBoundary <- st_difference(coastBuff,coast)  #Get boundary of analysis area
dataBoundary <- st_sfc(st_polygon(dataBoundary[[1]][[2]][1]),crs=st_crs(dataBoundary)) #Get rid of small islands


#Get file paths
files <-  data.frame(paths=dir(storage,recursive=TRUE,full.names = TRUE)) %>% 
  mutate(file=dir(storage,recursive=TRUE)) %>%
  filter(!grepl('combined',paths)) %>% 
  mutate(file=gsub('.*\\/','',file)) %>% 
  mutate(doy=as.numeric(gsub('\\.tif','',gsub('[AT]2014','',file)))) %>%
  dplyr::select(-file) %>% 
  mutate(platform=ifelse(grepl('terr',paths),'terra','aqua')) %>% 
  mutate(doy=ifelse(platform=='terra',doy-1,doy)) %>% #Terra starts 1 day earlier
  mutate(date=format(as.Date(paste0(doy,' 2014'),format='%j %Y'),format='%b %d')) %>%
  pivot_wider(names_from=platform,values_from = paths) %>%
  mutate(combined=gsub('\\/media.*\\/[AT]','',aqua)) %>% 
  mutate(combined=gsub('2014','2014_',combined)) %>% 
  mutate(combined=paste0(newfolder,'/',combined)) %>% 
  mutate(exists=file.exists(combined)) 

# #Fiddling around with dates
# as.Date(paste0(files$doy,'-2014'),format='%j-%Y')
# as.Date('91-2014',format='%j-%Y')
# format(as.Date(c('May 1 2014', 'Oct 1 2014'),format='%b %d %Y'),format='%j')

# #Create files of average between aqua and terra
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
#     print(paste0('Created file ',i,' of ',nrow(files)))
#   } else {
#     print(paste0('File ',i,' of ',nrow(files),' already exists'))
#   }
# }
# b <- Sys.time()
# b-a #Takes about 6 mins for 30 images

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
sDat2 <- cbind(sDat2,PCs) #Combine PCs with sDat2# 

#Add coordinate system
sDat2 <- sDat2 %>% st_as_sf(coords=c('x','y')) %>% st_set_crs(4326) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) #Louisiana offshore

Emean <- mean(unique(sDat2$E)); Nmean <- mean(unique(sDat2$N)) #Center values of E/N, used for scaling

sDat2 <- sDat2 %>% 
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% 
  st_transform(4326) %>% dplyr::select(-chlor_a:-sst) %>% unite(loc,E,N) %>% 
  mutate(loc=as.numeric(factor(loc)))

withinBuff <- sDat2 %>% st_intersects(.,dataBoundary) %>% sapply(.,function(x) length(x)>0) #Points in sDat2 that are outside of the buffer
sDat2 <- sDat2 %>% filter(withinBuff) #Filter out points outside of buffer
locLookup <- sDat2 %>% dplyr::select(loc:geometry) %>% unique() #Lookup table for locations
# ggplot(locLookup)+geom_sf()

rm(obs,sDat,withinBuff); gc()

# Fit PC models -----------------------------------

# #Approach using thin-plate splines

# library(parallel)
# cl <- makeCluster(15)
# #Takes about 10 minutes using 15 cores
# bFuns <- c('tp','tp'); kNum <- c(60,20); dNum <- c(2,1)
# a <- Sys.time()
# PCmod1 <- bam(PC1~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl) 
# PCmod2 <- bam(PC2~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod3 <- bam(PC3~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod4 <- bam(PC4~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod5 <- bam(PC5~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# PCmod6 <- bam(PC6~te(sN,sE,doy,bs=bFuns,k=kNum,d=dNum),data=sDat2,cluster=cl)
# par(mfrow=c(2,2)); gam.check(PCmod1); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
# par(mfrow=c(2,2)); gam.check(PCmod2); abline(0,1,col='red'); par(mfrow=c(1,1)) 
# par(mfrow=c(2,2)); gam.check(PCmod3); abline(0,1,col='red'); par(mfrow=c(1,1)) 
# par(mfrow=c(2,2)); gam.check(PCmod4); abline(0,1,col='red'); par(mfrow=c(1,1))
# par(mfrow=c(2,2)); gam.check(PCmod5); abline(0,1,col='red'); par(mfrow=c(1,1))
# par(mfrow=c(2,2)); gam.check(PCmod6); abline(0,1,col='red'); par(mfrow=c(1,1))
# Sys.time()-a
# stopCluster(cl)
# save(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6,file='./data/PCmods_mapping.RData')
# load('./data/PCmods_mapping.RData') #Load original tensor product smoothers

# #Approach using soap-film smoothers
# #Testing:
# 
# # #Buffer data area - knots can't be too close to boundary
# # dbTemp <- dataBoundary %>% st_transform(3401) %>% st_buffer(-30*1000) %>% st_transform(st_crs(dataBoundary))
# # #Create sampling zones for knots - idea: many close to the coast, fewer further away
# # zones <- lapply(c(75,150),function(d){
# #   coast[[1]][[1]] %>% st_polygon() %>% st_sfc(.,crs=st_crs(coast)) %>% #Only largest segment
# #     st_transform(3401) %>% st_buffer(d*1000) %>% st_sf() %>% st_cast('POLYGON') %>%
# #     mutate(area=as.numeric(st_area(.))) %>% filter(area==max(area)) %>%
# #     st_geometry() %>% st_transform(st_crs(dataBoundary)) %>%
# #     st_intersection(.,dbTemp)
# # })
# #
# # zones[[length(zones)+1]] <- st_difference(dbTemp,zones[[length(zones)+1-1]]) #Creates last zone
# # for(i in (length(zones)-1):2){
# #   zones[[i]] <- st_difference(zones[[i]],zones[[i-1]])
# # }
# # rm(dbTemp)
# #
# # #Looks OK
# #
# # # plot(knotLocs,add=TRUE,pch=19)
# # knotLocs <- mapply(function(z,N){
# #   z %>% st_transform(3401) %>% st_sample(size=N,type='hexagonal') %>% st_transform(st_crs(z))
# # },z=zones,N=c(30,20,20)) %>%
# #   do.call('c',.)
# #
# # plot(dataBoundary)
# # plot(zones[[1]],add=TRUE,col='green') #Looks OK
# # plot(zones[[2]],add=TRUE,col='red')
# # plot(zones[[3]],add=TRUE,col='blue')
# # plot(knotLocs,add=TRUE,pch=19)
# 
# # #Strips out knots that are too close - ideally this would use some other kind of approach like simulated annealing to maximize distance between points while holding them within sampling bands
# # minDist <- 30000 #20 km
# # removeKnot <- st_distance(knotLocs) %>% apply(.,1,function(x) min(x[x!=0])<minDist) #Should any knots be removed?
# # while(any(removeKnot)){
# #   knotLocs <- knotLocs[-which(removeKnot)[1],]
# #   removeKnot <- st_distance(knotLocs) %>% apply(.,1,function(x) min(x[x!=0])<minDist) #Checks if any points are closer than minDist
# # }
# 
# # # Alternative: constant knot locations
# # knotLocs <- dataBoundary %>% st_transform(3401) %>% st_buffer(-20*1000) %>%
# #   st_transform(st_crs(dataBoundary)) %>%
# #   st_transform(3401) %>% st_sample(size=100,type='hexagonal') %>% st_transform(st_crs(dataBoundary))
# # plot(dataBoundary)
# # plot(knotLocs,add=TRUE,pch=19)
# 
# #Alternative: knots along lines of buffer distances
# knotLocs <- mapply(function(d,N){
#   dbTemp <- dataBoundary %>% st_transform(3401) %>% st_buffer(-30*1000) %>% st_transform(st_crs(dataBoundary))
#     coast[[1]][[1]] %>% st_polygon() %>% st_sfc(.,crs=st_crs(coast)) %>% #Only largest segment
#     st_transform(3401) %>% st_buffer(d*1000) %>% st_sf() %>% st_cast('POLYGON') %>%
#     mutate(area=as.numeric(st_area(.))) %>% filter(area==max(area)) %>%
#     st_geometry() %>% st_transform(st_crs(dataBoundary)) %>%
#     st_cast('LINESTRING') %>% st_intersection(dbTemp,.) %>% #Transforms to line, gets intersection with dbTemp
#     st_cast('MULTILINESTRING') %>% st_line_merge() %>%  st_transform(3401) %>% #Converts to single line if needed
#     st_line_sample(n = N,type='regular') %>% #Samples N points along line
#     st_cast('POINT') %>%  st_transform(st_crs(dataBoundary)) #Transforms back to starting crs
# },d=c(50,125,225),N=c(30,20,10)) %>%  #Distances = 50,125,
#   do.call('c',.)
# 
# plot(dataBoundary) #Looks OK
# plot(knotLocs,add=TRUE,pch=19)
# length(knotLocs)
# 
# #Trying smaller example - single day
# tempDat <- sDat2 %>% filter(doy==211) %>% st_drop_geometry()
# 
# #Get knot and boundary locations on same scale
# 
# bound <- dataBoundary %>% st_transform(3401) %>%
#   st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   mutate(across(E:N,~round(.x))) %>% dplyr::select(sN,sE)
# bound <- list(list(sE=bound$sE,sN=bound$sN)) #Convert to list
# 
# kts <- knotLocs %>% st_transform(3401) %>%
#   st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   mutate(across(E:N,~round(.x))) %>% dplyr::select(sE,sN)
# 
# # with(bound[[1]],plot(sE,sN,type='l')) #Trying same thing with manually placed knots
# # kts <- locator() #Click on map for knots
# # names(kts) <- c('sE','sN')
# # points(kts$sE,kts$sN,pch=19,cex=0.8)
# # with(kts,inSide(bound,sE,sN)) #Check whether knots are inside boundary
# # in.out(bnd=bound,x=do.call('cbind',kts))
# 
# {a <- Sys.time() #Soap film
# PCmod1a <- gam(PC1~s(sE,sN,bs='so',xt=list(bnd=bound,nmax=100),k=30),
#               knots=kts,data=tempDat)
# Sys.time()-a}
# 
# {a <- Sys.time() #Thin-plate
# PCmod1b <- gam(PC1~s(sE,sN,bs='tp',k=60),data=tempDat)
# Sys.time()-a}
# 
# #Check results
# mae(PCmod1a); mae(PCmod1b)
# par(mfrow=c(3,1))
# plot(PCmod1a)
# plot(PCmod1b,scheme=2)
# pal <- colorRampPalette(c("blue", "red"))
# tempDat$order = findInterval(tempDat$PC1, sort(tempDat$PC1))
# with(tempDat,plot(sE,sN,col=pal(nrow(tempDat))[order],pch=19,cex=0.5))
# par(mfrow=c(1,2))
# plot(predict(PCmod1a),tempDat$PC1,main='Soap film'); abline(0,1,col='red')
# plot(predict(PCmod1b),tempDat$PC1,main='Thin-plate'); abline(0,1,col='red')
# par(mfrow=c(1,1))
# 
# #Model with tensor product smooth for days 
# 
# #Uses a subset of days: 201-221
# tempDat <- sDat2 %>% filter(doy>=201,doy<=221) %>% st_drop_geometry()
# 
# library(parallel)
# cl <- makeCluster(15)
# {a <- Sys.time()
# PCmod1a <- bam(PC1~te(sN,sE,doy,bs=c('sf','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,5),d=c(2,1))+
#                  te(sN,sE,doy,bs=c('sw','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,5),d=c(2,1)),
#               knots=kts,
#               data=tempDat,cluster=cl)
# Sys.time()-a} #Takes 22 seconds
# 
# {a <- Sys.time()
# PCmod1b <-bam(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,5),d=c(2,1)),
#               data=tempDat,cluster=cl)
# Sys.time()-a} #Takes 10 seconds
# stopCluster(cl)
# 
# #Check results
# mae(PCmod1a); mae(PCmod1b)
# par(mfrow=c(3,1))
# plot(PCmod1a,scheme=2)
# plot(PCmod1b,scheme=2)
# 
# par(mfrow=c(1,2))
# plot(predict(PCmod1a),tempDat$PC1,main='Soap film'); abline(0,1,col='red')
# plot(predict(PCmod1b),tempDat$PC1,main='Thin-plate'); abline(0,1,col='red')
# par(mfrow=c(1,1))
# 
# dataBoundary %>% st_transform(3401) %>%
#   st_sample(1000,type='hexagonal') %>%
#   st_sf() %>%
#   geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>%  #Louisiana offshore
#   # st_drop_geometry() %>%
#   st_transform(4326) %>%
#   expand_grid(.,doy=c(201,211,222)) %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   dplyr::select(-E,-N) %>%
#   mutate(predSoap=predict(PCmod1a,newdata=.),predTP=predict(PCmod1b,newdata=.)) %>%
#   pivot_longer(cols=c(predSoap,predTP)) %>%
#   ggplot()+geom_sf(aes(geometry=geometry,col=value))+
#   facet_grid(doy~name)


# #Full model fits
# # Knots along lines of buffer distances
# knotLocs <- mapply(function(d,N){
#   dbTemp <- dataBoundary %>% st_transform(3401) %>% st_buffer(-30*1000) %>% st_transform(st_crs(dataBoundary))
#   coast[[1]][[1]] %>% st_polygon() %>% st_sfc(.,crs=st_crs(coast)) %>% #Only largest segment
#     st_transform(3401) %>% st_buffer(d*1000) %>% st_sf() %>% st_cast('POLYGON') %>%
#     mutate(area=as.numeric(st_area(.))) %>% filter(area==max(area)) %>%
#     st_geometry() %>% st_transform(st_crs(dataBoundary)) %>%
#     st_cast('LINESTRING') %>% st_intersection(dbTemp,.) %>% #Transforms to line, gets intersection with dbTemp
#     st_cast('MULTILINESTRING') %>% st_line_merge() %>%  st_transform(3401) %>% #Converts to single line if needed
#     st_line_sample(n = N,type='regular') %>% #Samples N points along line
#     st_cast('POINT') %>%  st_transform(st_crs(dataBoundary)) #Transforms back to starting crs
# },
# d=c(50,150,225), #Distances 
# N=c(25,15,8) #Points within each distance
# ) %>%  
#   do.call('c',.)
# 
# #Looks OK
# plot(dataBoundary)
# plot(knotLocs,add=TRUE,pch=19,cex=0.5)
# length(knotLocs)
# 
# #Get knot and boundary locations on same scale
# bound <- dataBoundary %>% st_transform(3401) %>%
#   st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   mutate(across(E:N,~round(.x))) %>% dplyr::select(sN,sE)
# bound <- list(list(sE=bound$sE,sN=bound$sN)) #Convert to list
# 
# kts <- knotLocs %>% st_transform(3401) %>%
#   st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   mutate(across(E:N,~round(.x))) %>% dplyr::select(sE,sN)
# 
# library(parallel)
# cl <- makeCluster(15)
# 
# #Per model: ~10 mins for 60 knots, 30 boundary loops, 15 layers for tensor product
# modForm <- "999~te(sN,sE,doy,bs=c('sf','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,20),d=c(2,1))+
#   te(sN,sE,doy,bs=c('sw','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,20),d=c(2,1))"
# a <- Sys.time()
# PCmod1 <- bam(formula(gsub('999','PC1',modForm)),knots=kts,data=sDat2,cluster=cl)
# Sys.time()-a; beep()
# PCmod2 <- bam(formula(gsub('999','PC2',modForm)),knots=kts,data=sDat2,cluster=cl)
# PCmod3 <- bam(formula(gsub('999','PC3',modForm)),knots=kts,data=sDat2,cluster=cl)
# PCmod4 <- bam(formula(gsub('999','PC4',modForm)),knots=kts,data=sDat2,cluster=cl)
# PCmod5 <- bam(formula(gsub('999','PC5',modForm)),knots=kts,data=sDat2,cluster=cl)
# PCmod6 <- bam(formula(gsub('999','PC6',modForm)),knots=kts,data=sDat2,cluster=cl)
# stopCluster(cl)
# # # Test against regular thin-plate spline
# # {a <- Sys.time() #5 mins
# #   PCmod1b <-bam(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(60,20),d=c(2,1)),
# #                 data=sDat2,cluster=cl)
# #   Sys.time()-a}
# # stopCluster(cl)
# # # mae(PCmod1); mae(PCmod1b) #Soap film smoother slightly better
# save(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6,file='./data/PCmodsSoap_mapping.RData')

load('./data/PCmodsSoap_mapping.RData') #Load soap film smoothers

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


#Goal: min/mean predictions for all points on following day segments

dateRanges <- read.csv('./data/mappingDateRanges.csv') %>% 
  mutate(across(c(startDay,endDay), ~as.numeric(format(as.Date(paste0(.x,' 2014'),format='%B %d %Y'),format='%j')),.names='{.col}_doy'))
dateLabs <- with(dateRanges,paste0(startDay,' - ',endDay)) #Labels for date ranges
dateRanges <- with(dateRanges,data.frame(slice=rep(slice,nDays),doy=min(startDay_doy):max(endDay_doy)))

d <- dateRanges$doy #All days
l <- sort(unique(locLookup$loc)) #All unique locations

newDF <- expand_grid(doy=d,loc=l) %>% #Days/locations to predict at
  mutate(doy2=doy,doy=doy-daylag) %>% #Actual day = doy2, lagged day to get data from = doy
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locLookup),by='loc')

#Separates into 2 df, one with values and one with NAs (quicker predictions)
newDF1 <- newDF %>% filter(!is.na(PC1))  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% dplyr::select(-contains('PC')) #DF with NAs
cl <- makeCluster(15)
a <- Sys.time() #Takes about 4 mins
newDF2 <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){ #Make predictions in NA spaces
  require(mgcv); predict.gam(x,newdata=N)
  },N=newDF2) %>% 
  set_names(paste0('PC',1:6)) %>% bind_cols() %>% bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN)
Sys.time()-a
stopCluster(cl)
newDF <- bind_rows(newDF1,newDF2) %>% arrange(loc,doy) %>% mutate(predDO=predict(m1,.)) %>% dplyr::select(-sE,-sN) 
rm(newDF1,newDF2) #Cleanup

#Data for maps
mapDat <- newDF %>% 
  left_join(.,dateRanges,by='doy') %>% 
  group_by(slice,loc) %>% summarize(minDO=min(predDO),meanDO=mean(predDO)) %>% 
  ungroup() %>% 
  mutate(across(c(minDO,meanDO),~cut(.x,breaks=c(min(.x,na.rm=TRUE),1,2,3,4,5,max(.x,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE))) %>% 
  mutate(slice=factor(slice,lab=dateLabs)) %>% 
  left_join(.,locLookup,by='loc') 

p1 <- ggplot(mapDat)+
  geom_sf(aes(col=minDO,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Lagged-linear - min(DO)',col='DO (mg/L)')
p2 <- ggplot(mapDat)+
  geom_sf(aes(col=meanDO,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Lagged-linear - mean(DO)',col='DO (mg/L)')

ggarrange(p1,p2,ncol=1)

  
#Get predictions from FR model ------------------------------



load('./data/funRegMod.Rdata')



#Possible prediction doy range | lag = days 182 - 243 (Jul 1 - Aug 31)

#Goal: 
# July 1-7, 8-14, 15-21, 21-28, 29-Aug 4, 5-11, 12-18,  19-24, 24-30.

# "July 01" "July 08" "July 15" "July 22" "July 30"
d <- c(182,189,196,203,211,243)
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
  set_names(paste0('PC',1:6)) %>% bind_cols() %>% 
  bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN)
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
