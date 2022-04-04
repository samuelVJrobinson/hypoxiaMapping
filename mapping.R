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
library(missMDA)
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

{a <- Sys.time()
  sDat <- lapply(files$combined,function(datPath){
    cNames <- names(brick(files$aqua[1])) #Channel names
    sD <- brick(datPath) #Read data
    names(sD) <- cNames #Add names
    sD2 <- rasterToPoints(sD) #Convert to point dataframe
    nMissing <- apply(sD2[,!grepl('(x|y)',colnames(sD2))],1,function(x) sum(is.na(x))) #Proportion of missing values in each row (cell)
    # sD3 <- sD2[nMissing<5,] #Retains records with 4 or more data values
    # if(nrow(sD3)==0) sD3 <- rbind(sD3,rep(NA,ncol(sD3)))
    sD4 <- data.frame(doy=strsplit(x = datPath,split = c("(\\_|\\.)"))[[1]][2], #Day of year
                      sD2,numNA=nMissing) 
    return(sD4)})
Sys.time()-a} #Takes about 2 mins
names(sDat) <- sapply(sDat,function(x) x$doy[1])

# #Preview day 2
# sDat[[2]] %>% #sst on day 2
#   st_as_sf(coords=c('x','y')) %>% st_set_crs(4326) %>%
#   ggplot()+
#   geom_sf(aes(col=sst))+
#   geom_sf(data=dataBoundary,fill=NA)

sDat2 <- do.call('rbind',sDat) #Combine into single DF 

# lwrLimits <- apply(dplyr::select(sDat2,chlor_a:sst),2,function(x) min(x[x>0],na.rm=TRUE)) #Minimum nonzero values
lwrLimits <- apply(dplyr::select(sDat2,chlor_a:sst),2,function(x) min(x,na.rm=TRUE)) #Minimum values
names(lwrLimits) <- c('chlor_a','nflh','poc','Rrs_412','Rrs_443','Rrs_469','Rrs_488','Rrs_531','Rrs_547','Rrs_555','Rrs_645','Rrs_667','Rrs_678','sst')

sDat2 <- sDat2 %>% 
  # # #Rescales negative values be above 0.95*minimum positive value
  # mutate(across(chlor_a:sst,~ifelse(.x<lwrLimits[names(lwrLimits)==cur_column()],lwrLimits[names(lwrLimits)==cur_column()]*0.95,.x)))
  # Adds lowest value + 10% if lwrLimits < 0
  mutate(across(names(lwrLimits[lwrLimits<0]),~.x-(min(.x,na.rm=TRUE)*1.1))) %>% 
  mutate(across(chlor_a:sst,log)) #Log-scale variables

sDatMat <- sDat2 %>% dplyr::select(chlor_a:sst) %>% as.matrix() #Data in matrix form
sDat2$numNA <- apply(sDatMat,1,function(x) sum(is.na(x))) #Get number of NA values
round(table(sDat2$numNA)/nrow(sDat2),4) #~50% cells missing all data, ~40% complete
locs <- dplyr::select(sDat2,x,y) %>% unique() %>% mutate(locID=1:n()) #Locations only
datMissing <- sDat2$numNA==13 #Locations missing all data

# sDat2Miss <- sDat2 %>% filter(datMissing) #Completely missing 
sDat2SomeMiss <- sDat2 %>% filter(numNA!=13,numNA!=0) #Some non missing
sDat2NoMiss <- sDat2 %>% filter(numNA==0) #No missing

sDatMat <- sDat2SomeMiss %>% dplyr::select(chlor_a:sst) %>% as.matrix() #Partly missing data in matrix form

# n_components <- estim_ncpPCA(sDatMat, verbose = TRUE,method="Regularized",method.cv="gcv",ncp.min=1,ncp.max=13) #Takes about 5 mins
# plot(1:13,n_components$criterion,xlab='Number of Dimensions',ylab='GCV Criterion',type='b',pch=19) #Looks like about 7 components is OK for prediction

#Same as imputePCA, but multicore. Splits X into slices of nPer each, passes to cluster cl, reassembles
imputePCA_multi <- function(X,ncp=2,scale=TRUE,method='Regularized',nPer=10000,ncore=10){
  # X <- sDatMat #Debugging
  # nPer <- 10000
  # ncore <- 10
  # ncp <- 3
  # scale <- TRUE
  # method <- 'Regularized'
  
  #Assign rows to sample strata
  sampNum <- c(rep(1:(nrow(X) %/% nPer),each=nPer), #Regular samples
    rep((nrow(X) %/% nPer + 1),(nrow(X) %% nPer))) #"Leftovers"
  
  sampOrd <- sample(1:length(sampNum))
  sampNum <- sampNum[sampOrd] #Randomize
  
  rNames <- rownames(X) #Save rownames
  rownames(X) <- 1:nrow(X)
  
  IP <- function(x,ncp,scale,method) imputePCA(x,ncp,scale,method)$completeObs #Convenience function to pull only 1 matrix from output
  
  print('Breaking into separate matrices')
  X <- lapply(1:max(sampNum),function(x) X[sampNum==x,]) #Break into separate matrices
  if(ncore>1){
    library(parallel)
    cl <- makeCluster(ncore)
    print(paste('Running imputePCA across',ncore,'clusters'))
    clusterEvalQ(cl,library(missMDA)) #Loads missMDA on each cluster
    for(i in 1:(ceiling(length(X)/ncore))){ #Break into sets. Tends to hang on large 
      use <- ((i-1)*ncore+1):pmin(length(X),(i*ncore)) #Which sets of matrices to use
      X[use] <- parLapply(cl=cl,X=X[use],fun=IP,ncp=ncp,scale=scale,method=method) #Run imputePCA on separate matrices
      cat('.') #Something to look at
      # X <- parLapply(cl=cl,X=X,fun=IP,ncp=ncp,scale=scale,method=method) #Run imputePCA on separate matrices
    }
    stopCluster(cl); gc()
  } else {
    # X <- lapply(X,IP,ncp=ncp,scale=scale,method=method) #Run imputePCA on separate matrices
    times <- rep(NA,length(X))
    for(i in 1:length(X)){
      a <- Sys.time()
      X[[i]] <- IP(X[[i]],ncp=ncp,scale=scale,method=method)
      b <- Sys.time()
      times[i] <- difftime(b,a,units='secs')
      print(paste0('Finished ',i,'. Time: ',round(times[i],2),'. Average time: ',round(mean(times,na.rm=TRUE),2), 
                   '. Estimated remaining (mins): ',round(mean(times,na.rm=TRUE),2)*(length(X)-i)/60 ))
    } 
  }
  X <- do.call(rbind,X) #Recombine into single matrix
  
  if(length(rNames)!=nrow(na.omit(X))) stop('NA values still present in dataframe')
  
  X <- X[order(as.numeric(rownames(X))),] #Reorder
  rownames(X) <- rNames #Reset rownames
  
  print('Done')
  return(X)
}

#~3 mins
sDatMat_imputed <- imputePCA_multi(sDatMat,ncp=7,scale=TRUE,method='Regularized',nPer = 10000,ncore = 15) #Impute missing data using 7 dimensions

sDat2SomeMiss[which(names(sDat2SomeMiss)=='chlor_a'):which(names(sDat2SomeMiss)=='sst')] <- sDatMat_imputed #Rejoin

rm(sDatMat_imputed); gc() #Cleanup

#Put back into a single dataframe again
sDat2 <- bind_rows(sDat2SomeMiss,sDat2NoMiss) %>% dplyr::select(-numNA) %>% arrange(doy,x,y)
rm(sDat2SomeMiss,sDat2NoMiss,sDatMat,datMissing); gc()

sDatMat <- sDat2 %>% dplyr::select(chlor_a:sst) %>% as.matrix() #Data in matrix form again

useThese <- unname(which(apply(sDatMat,1,function(x) !any(is.na(x))))) #Rows without missing values

(pca1 <- prcomp(sDatMat[useThese,],center = TRUE, scale. = TRUE)) #Principle components
pca1$center
pca1$scale

#Looks like 5 PCs get 98.6% of var
data.frame(pc=1:length(pca1$sdev),cVar=cumsum(pca1$sdev^2)/sum(pca1$sdev^2)) %>%
    ggplot(aes(x=pc,y=cVar))+geom_point()+geom_line()+
    labs(x='Principle Component',y='Cumulative Variance')+
    geom_hline(yintercept = 0.95,col='red',linetype='dashed')+
  scale_x_continuous(n.breaks=length(pca1$sdev),minor_breaks = NULL)

#Factor loadings
pca1$rotation[,1:5] %>% data.frame() %>% rownames_to_column(var='var') %>%
  pivot_longer(PC1:PC5) %>%
  mutate(name=factor(name,labels=paste0('PC',1:5,': ',round(pca1$sdev[1:5]^2/sum(pca1$sdev^2),3)*100,'% Variance'))) %>%
  mutate(var=factor(var,levels=rev(unique(var)))) %>%
  ggplot()+ geom_col(aes(y=var,x=value))+geom_vline(xintercept = 0,col='red',linetype='dashed')+
  facet_wrap(~name)+
  labs(x='Loading',y=NULL,title='Factor loadings for Principle Components 1-5')

#Matrix to store PC values
pcaMat <- matrix(NA,nrow = nrow(sDatMat), ncol=5,dimnames = list(rownames(sDatMat),paste0('PC',1:5)))

pcaMat[useThese,] <- pca1$x[,1:ncol(pcaMat)] #Store PC values in matrix

sDat2 <- cbind(sDat2,pcaMat) #Combine PCs with sDat2 

rm(sDatMat,pcaMat); gc()

#Add coordinate system - takes a minute or so
sDat2 <- sDat2 %>% st_as_sf(coords=c('x','y')) %>% 
  dplyr::select(-chlor_a:-sst) %>% #Remove raw data, keep PCs
  st_set_crs(4326) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) #Louisiana offshore

Emean <- mean(unique(sDat2$E)); Nmean <- mean(unique(sDat2$N)) #Center values of E/N, used for scaling

sDat2 <- sDat2 %>% 
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% 
  st_transform(4326) %>% unite(loc,E,N) %>% 
  mutate(loc=as.numeric(factor(loc)))

withinBuff <- sDat2 %>% st_intersects(.,dataBoundary) %>% sapply(.,function(x) length(x)>0) #Points in sDat2 that are outside of the buffer
sDat2 <- sDat2 %>% filter(withinBuff) #Filter out points outside of buffer

rm(withinBuff,useThese); gc()

# PC1-3 changes through the year
# NOTE: Many locations/times are completely missing
sDat2 %>% 
  filter(doy=='071'|doy=='141'|doy=='211') %>% 
  pivot_longer(cols=c(PC1:PC3)) %>% 
  ggplot()+
  geom_sf(data=dataBoundary,fill=NA)+
  geom_sf(aes(col=value))+
  facet_grid(doy~name)

# ggplot()+geom_sf(dat=dataBoundary)+geom_sf(dat=locLookup)

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

#Approach using soap-film smoothers
#Testing:
# #Buffer data area - knots can't be too close to boundary
# dbTemp <- dataBoundary %>% st_transform(3401) %>% st_buffer(-30*1000) %>% st_transform(st_crs(dataBoundary))
# #Create sampling zones for knots - idea: many close to the coast, fewer further away
# zones <- lapply(c(75,150),function(d){
#   coast[[1]][[1]] %>% st_polygon() %>% st_sfc(.,crs=st_crs(coast)) %>% #Only largest segment
#     st_transform(3401) %>% st_buffer(d*1000) %>% st_sf() %>% st_cast('POLYGON') %>%
#     mutate(area=as.numeric(st_area(.))) %>% filter(area==max(area)) %>%
#     st_geometry() %>% st_transform(st_crs(dataBoundary)) %>%
#     st_intersection(.,dbTemp)
# })
#
# zones[[length(zones)+1]] <- st_difference(dbTemp,zones[[length(zones)+1-1]]) #Creates last zone
# for(i in (length(zones)-1):2){
#   zones[[i]] <- st_difference(zones[[i]],zones[[i-1]])
# }
# rm(dbTemp)
#
# #Looks OK
#
# # plot(knotLocs,add=TRUE,pch=19)
# knotLocs <- mapply(function(z,N){
#   z %>% st_transform(3401) %>% st_sample(size=N,type='hexagonal') %>% st_transform(st_crs(z))
# },z=zones,N=c(30,20,20)) %>%
#   do.call('c',.)
#
# plot(dataBoundary)
# plot(zones[[1]],add=TRUE,col='green') #Looks OK
# plot(zones[[2]],add=TRUE,col='red')
# plot(zones[[3]],add=TRUE,col='blue')
# plot(knotLocs,add=TRUE,pch=19)

# #Strips out knots that are too close - ideally this would use some other kind of approach like simulated annealing to maximize distance between points while holding them within sampling bands
# minDist <- 30000 #20 km
# removeKnot <- st_distance(knotLocs) %>% apply(.,1,function(x) min(x[x!=0])<minDist) #Should any knots be removed?
# while(any(removeKnot)){
#   knotLocs <- knotLocs[-which(removeKnot)[1],]
#   removeKnot <- st_distance(knotLocs) %>% apply(.,1,function(x) min(x[x!=0])<minDist) #Checks if any points are closer than minDist
# }

# # Alternative: constant knot locations
# knotLocs <- dataBoundary %>% st_transform(3401) %>% st_buffer(-20*1000) %>%
#   st_transform(st_crs(dataBoundary)) %>%
#   st_transform(3401) %>% st_sample(size=100,type='hexagonal') %>% st_transform(st_crs(dataBoundary))
# plot(dataBoundary)
# plot(knotLocs,add=TRUE,pch=19)

# #Alternative: knots along lines of buffer distances
# 
# #Small example - single day
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
# ggplot() + #Looks OK
#   geom_sf(data=dataBoundary,fill=NA) +
#   geom_sf(data=slice_sample(sDat2,n=5000),alpha=0.3)+
#   geom_sf(data=knotLocs,col='red') 
# 
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
# 
# #Small example - tensor product smooth across days
# #Uses a subset of days: 201-221
# tempDat <- sDat2 %>% filter(doy>=201,doy<=221) %>% st_drop_geometry()
# tempGeom <- sDat2 %>% filter(doy>=201,doy<=221) %>% st_geometry()
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
# sampRows <- sample(1:nrow(tempDat),size = 5000)
# plot(predict(PCmod1a)[sampRows],tempDat$PC1[sampRows],main='Soap film'); abline(0,1,col='red')
# plot(predict(PCmod1b)[sampRows],tempDat$PC1[sampRows],main='Thin-plate'); abline(0,1,col='red')
# par(mfrow=c(1,1))
# 
# #Predicted values
# dataBoundary %>% st_transform(3401) %>%
#   st_sample(1000,type='hexagonal') %>%
#   st_sf() %>%
#   geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>%  #Louisiana offshore
#   # st_drop_geometry() %>%
#   st_transform(4326) %>%
#   expand_grid(.,doy=c(201,211,221)) %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   dplyr::select(-E,-N) %>%
#   mutate(predSoap=predict(PCmod1a,newdata=.),predTP=predict(PCmod1b,newdata=.)) %>%
#   pivot_longer(cols=c(predSoap,predTP)) %>%
#   ggplot()+
#   geom_sf(aes(geometry=geometry,col=value))+
#   facet_grid(doy~name)
# 
# #Error
# tempDat %>% 
#   mutate(errSoap=predict(PCmod1a)-tempDat$PC1,errTP=predict(PCmod1b)-tempDat$PC1) %>% 
#   st_sf(geom=tempGeom) %>% 
#   filter(doy==201|doy==211|doy==221) %>% 
#   # dplyr::select(doy,contains('err')) %>% 
#   pivot_longer(cols=c(errSoap,errTP)) %>%
#   ggplot()+
#   geom_sf(aes(geometry=geom,col=value))+
#   facet_grid(doy~name)+
#   scale_color_distiller(palette = 'Spectral')


#Full model fits
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

ggplot() + #Looks OK
  geom_sf(data=dataBoundary,fill=NA) +
  geom_sf(data=slice_sample(sDat2,n=5000),alpha=0.3)+
  geom_sf(data=knotLocs,col='red')

#Get knot and boundary locations on same scale
bound <- dataBoundary %>% st_transform(3401) %>%
  st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% dplyr::select(sN,sE)
bound <- list(list(sE=bound$sE,sN=bound$sN)) #Convert to list

kts <- knotLocs %>% st_transform(3401) %>%
  st_coordinates() %>% data.frame() %>% rename(E=X,N=Y) %>%
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% dplyr::select(sE,sN)

#Goal: split up sDat2 into manageable chunks for soap film smoothing

sDat2 %>% st_drop_geometry() %>% group_by(doy) %>% 
  summarize(#propNA=sum(is.na(PC1))/n(),
            n=n())#,nData=sum(!is.na(PC1))) %>% 
  mutate(doy=as.numeric(doy)) %>% 
  ggplot(aes(x=doy,y=propNA))+geom_point()
  


library(parallel)
cl <- makeCluster(15)

#Per model: ~10 mins for 60 knots, 30 boundary loops, 20 layers for tensor product 
modForm <- "999~te(sN,sE,doy,bs=c('sf','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,30),d=c(2,1))+
  te(sN,sE,doy,bs=c('sw','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(30,30),d=c(2,1))"
a <- Sys.time()
PCmod1 <- bam(formula(gsub('999','PC1',modForm)),knots=kts,data=sDat2,cluster=cl)
Sys.time()-a; beep()
PCmod2 <- bam(formula(gsub('999','PC2',modForm)),knots=kts,data=sDat2,cluster=cl)
PCmod3 <- bam(formula(gsub('999','PC3',modForm)),knots=kts,data=sDat2,cluster=cl)
PCmod4 <- bam(formula(gsub('999','PC4',modForm)),knots=kts,data=sDat2,cluster=cl)
PCmod5 <- bam(formula(gsub('999','PC5',modForm)),knots=kts,data=sDat2,cluster=cl)
PCmod6 <- bam(formula(gsub('999','PC6',modForm)),knots=kts,data=sDat2,cluster=cl)
stopCluster(cl)
# # Test against regular thin-plate spline
# {a <- Sys.time() #5 mins
#   PCmod1b <-bam(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(60,20),d=c(2,1)),
#                 data=sDat2,cluster=cl)
#   Sys.time()-a}
# stopCluster(cl)
# # mae(PCmod1); mae(PCmod1b) #Soap film smoother slightly better
save(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6,file='./data/PCmodsSoap_mapping.RData')

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


#Get min/mean predictions for all points on following day segments

dateRanges <- read.csv('./data/mappingDateRanges.csv') %>% 
  mutate(across(c(startDay,endDay), ~as.numeric(format(as.Date(paste0(.x,' 2014'),format='%B %d %Y'),format='%j')),.names='{.col}_doy'))
dateLabs <- with(dateRanges,paste0(startDay,' - ',endDay)) #Labels for date ranges
dateRanges <- with(dateRanges,data.frame(slice=rep(slice,nDays),doy=min(startDay_doy):max(endDay_doy)))

d <- dateRanges$doy #All days
l <- sort(unique(locLookup$loc)) #All unique locations
daylag <- 12 #Use 12-day lag

newDF <- expand_grid(doy=d,loc=l) %>% #Days/locations to predict at
  mutate(doy2=doy,doy=doy-daylag) %>% #doy = lagged day to get DO predictions from, doy 2 = Actual day 
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locLookup),by='loc')

#Separates into 2 df, one with values and one with NAs (quicker predictions)
newDF1 <- newDF %>% filter(!is.na(PC1))%>% mutate(imputed=FALSE)  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% dplyr::select(-contains('PC')) #DF with NAs
cl <- makeCluster(15)
{a <- Sys.time() #Takes about 4 mins
newDF2 <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){ #Make predictions in NA spaces
  require(mgcv); predict.gam(x,newdata=N)
  },N=newDF2) %>% 
  set_names(paste0('PC',1:6)) %>% bind_cols() %>% bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN) %>% mutate(imputed=TRUE)
Sys.time()-a}
stopCluster(cl)
newDF <- bind_rows(newDF1,newDF2) %>% 
  arrange(loc,doy) %>% mutate(predDO=predict(m1,.)) %>% 
  dplyr::select(-sE,-sN,-doy) %>% rename('doy'='doy2')
# rm(newDF1,newDF2) #Cleanup

#Data for maps
mapDat <- newDF %>% 
  left_join(.,dateRanges,by='doy') %>% 
  group_by(slice,loc) %>% summarize(minDO=min(predDO),meanDO=mean(predDO),propImputed=sum(imputed)/n()) %>% 
  ungroup() %>% 
  mutate(across(c(minDO,meanDO),~cut(.x,breaks=c(min(.x,na.rm=TRUE),1,2,3,4,5,max(.x,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE),
                .names='{.col}_factor')) %>% 
  mutate(slice=factor(slice,lab=dateLabs)) %>% 
  left_join(locLookup,.,by='loc') %>% dplyr::select(-loc,-sE,-sN)

#Save a copy for YingJie
LLmapDat <- mapDat %>% 
  relocate(slice,minDO,meanDO,minDO_factor,meanDO_factor,propImputed) %>% 
  geom2cols()
save(LLmapDat,file = './data/LLmapDat.Rdata')

#Full copy for YingJie
LLmapDatFull <- newDF %>% left_join(.,locLookup,by='loc') %>% 
  dplyr::select(doy,imputed,predDO,geometry) %>% st_as_sf() %>% 
  geom2cols()
save(LLmapDatFull,file = './data/LLmapDatFull.Rdata')

#Make figure
p1 <- ggplot(mapDat)+
  geom_sf(aes(col=minDO_factor,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Lagged-linear - min(DO)',col='DO (mg/L)')
p2 <- ggplot(mapDat)+
  geom_sf(aes(col=meanDO_factor,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Lagged-linear - mean(DO)',col='DO (mg/L)')
p <- ggarrange(p1,p2,ncol=1)

ggsave(p,filename = './figures/mapLLpred.png',width=12,height=15)

# Get predictions from FR model ------------------------------

load('./data/funRegMod.Rdata')

dateRanges <- read.csv('./data/mappingDateRanges.csv') %>% #Get dates to predict on
  mutate(across(c(startDay,endDay), ~as.numeric(format(as.Date(paste0(.x,' 2014'),format='%B %d %Y'),format='%j')),.names='{.col}_doy'))
dateLabs <- with(dateRanges,paste0(startDay,' - ',endDay)) #Labels for date ranges
dateRanges <- with(dateRanges,data.frame(slice=rep(slice,nDays),doy=min(startDay_doy):max(endDay_doy)))

d <- dateRanges$doy #All days
l <- sort(unique(locLookup$loc)) #All unique locations
daylag <- 30 #30-day lag

#Dataframe of values to choose from
newDF <- expand_grid(doy=unique(do.call('c',lapply(d,function(x) x-(daylag:0)))),loc=l) %>% 
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locLookup),by='loc')
newDF1 <- newDF %>% filter(!is.na(PC1)) %>% mutate(imputed=FALSE)  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% dplyr::select(-contains('PC')) #DF with NAs
cl <- makeCluster(15)
{a <- Sys.time() #Takes about 4.6 mins for all 2.4 mil points
newDF2 <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){require(mgcv); predict.gam(x,newdata=N)},N=newDF2) %>% 
  set_names(paste0('PC',1:6)) %>% bind_cols() %>% 
  bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN) %>% mutate(imputed=TRUE)
Sys.time()-a}; beep()
stopCluster(cl)
newDF <- bind_rows(newDF1,newDF2) %>% arrange(loc,doy)
rm(newDF1,newDF2) #Cleanup

#Make list of matrices
datList <- with(expand.grid(l=l,d=d),list(loc=l,doy=d))
datList$dayMat <- outer(rep(1,length(datList$doy)),0:daylag)

#For each day, create dataframe of all locations and all lagged PC values - not super efficient but it works
nd <- lapply(1:length(d),function(i){ #Takes about a minute
  chooseThese <- (newDF$doy %in% (d[i]-daylag):d[i]) & (newDF$loc %in% l)
  newDF[chooseThese,] %>% data.frame()
}) %>% bind_rows() 

pcaMats <- lapply(which(grepl('PC',names(nd))),function(j){
  return(matrix(nd[,j],ncol=ncol(datList$dayMat),byrow=TRUE)[,ncol(datList$dayMat):1])
}) %>% set_names(paste0('pcaMat',1:6))

datList <- c(datList,pcaMats)
datList$isImputed <- matrix(nd[,which(grepl('imputed',names(nd)))],ncol=ncol(datList$dayMat),byrow=TRUE) #Which values are imputed?

{a <- Sys.time() #Takes 4.6 mins to generate predictions for each day/location
  mapDat <- with(datList,data.frame(loc=loc,doy=doy,predDO=predict(bWatMod,newdata=datList))) 
  Sys.time()-a}

#Gets proportion of imputed data used in each prediction
mapDat$propImputed <- apply(datList$isImputed,1,mean)

#Full copy for YingJie
FRmapDatFull <- mapDat %>% left_join(.,locLookup,by='loc') %>% 
  dplyr::select(doy,propImputed,predDO,geometry) %>% st_as_sf() %>% 
  geom2cols()
save(FRmapDatFull,file = './data/FRmapDatFull.Rdata')

mapDat <- mapDat %>% 
  left_join(.,dateRanges,by='doy') %>% 
  group_by(slice,loc) %>% summarize(minDO=min(predDO),meanDO=mean(predDO),propImputed=mean(propImputed)) %>% 
  ungroup() %>% 
  mutate(across(c(minDO,meanDO),~cut(.x,breaks=c(min(.x,na.rm=TRUE),1,2,3,4,5,max(.x,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE),
                .names='{.col}_factor')) %>% 
  mutate(slice=factor(slice,lab=dateLabs)) %>% 
  left_join(locLookup,.,by='loc') %>% dplyr::select(-loc,-sE,-sN)

#Save data for Yingie
FRmapDat <- mapDat %>% 
  relocate(slice,minDO,meanDO,minDO_factor,meanDO_factor,propImputed) %>% 
  mutate(propImputed=NA) %>% 
  geom2cols()
save(FRmapDat,file = './data/FRmapDat.Rdata')

#Make figure
p1 <- ggplot(mapDat)+
  geom_sf(aes(col=minDO_factor,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Functional Regression - min(DO)',col='DO (mg/L)')
p2 <- ggplot(mapDat)+
  geom_sf(aes(col=meanDO_factor,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Functional Regression - mean(DO)',col='DO (mg/L)')
p <- ggarrange(p1,p2,ncol=1)
ggsave(p,filename = './figures/mapFRpred.png',width=12,height=15)

# Convert predictions from dataframe to raster -----------------------------

library(sp)
library(raster)

#Get info needed for transforming points to rasters
tempRaster <- brick("./data/2014_121.tif") #Read in first combined tif file
tempRes <- res(tempRaster) #Resolution
tempCRS <- crs(tempRaster) #CRS
tempExtent <- extent(tempRaster) #Extent

#Load in predictions dataframe
load('./data/FRmapDat.Rdata') #This isn't located on the Git, so download from Google Drive first

#Not sure what your workflow is, but say you only wanted a raster for the first slice:
tempDF <- filter(FRmapDat,slice=='June 1 - June 10') #Choose only first slice
tempCoords <- tempDF[,c('lon','lat')]
tempDF <- tempDF %>% dplyr::select(minDO,meanDO) #Get rid of extra info

tempDF <- SpatialPointsDataFrame(coords=tempCoords,data=tempDF,proj4string=tempCRS) #Convert to sp object

exampleRaster <- raster(crs = tempCRS, vals = NA, resolution = tempRes, ext = tempExtent) %>%
  rasterize(tempDF, .)

plot(exampleRaster) #Works

#Now that it's been converted, you can save it in whatever form you want using writeRaster()


