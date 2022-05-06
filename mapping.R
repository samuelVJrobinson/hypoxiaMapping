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

load('./data/all2014_3.Rdata')
rm(wDat,sDat,surfWDat)#,locIndex,bottomWDat)

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

# #Get file paths
# files <-  data.frame(paths=dir(storage,recursive=TRUE,full.names = TRUE)) %>%
#   mutate(file=dir(storage,recursive=TRUE)) %>%
#   filter(!grepl('combined',paths)) %>%
#   mutate(file=gsub('.*\\/','',file)) %>%
#   mutate(doy=as.numeric(gsub('\\.tif','',gsub('[AT]2014','',file)))) %>%
#   dplyr::select(-file) %>%
#   mutate(platform=ifelse(grepl('terr',paths),'terra','aqua')) %>%
#   mutate(doy=ifelse(platform=='terra',doy-1,doy)) %>% #Terra starts 1 day earlier
#   mutate(date=format(as.Date(paste0(doy,' 2014'),format='%j %Y'),format='%b %d')) %>%
#   filter(doy>65) %>% #Only use March 6th and beyond
#   pivot_wider(names_from=platform,values_from = paths) %>%
#   mutate(combined=gsub('\\/media.*\\/[AT]','',aqua)) %>%
#   mutate(combined=gsub('2014','2014_',combined)) %>%
#   mutate(combined=paste0(newfolder,'/',combined)) %>%
#   mutate(exists=file.exists(combined))
# 
# # #Fiddling around with dates
# # as.Date(paste0(files$doy,'-2014'),format='%j-%Y')
# # as.Date('91-2014',format='%j-%Y')
# # format(as.Date(c('May 1 2014', 'Oct 1 2014'),format='%b %d %Y'),format='%j')
# 
# # #Create files of average between aqua and terra
# # a <- Sys.time()
# # for(i in 1:nrow(files)){
# #   if(!files$exists[i]){
# #     aquaDat <- brick(files$aqua[i])
# #     terraDat <- brick(files$terra[i])
# #
# #     newDat <- lapply(1:dim(aquaDat)[3],function(j) calc(stack(aquaDat[[j]],terraDat[[j]]),fun=mean,na.rm=TRUE))
# #     newDat <- brick(newDat)
# #     names(newDat) <- names(aquaDat)
# #     writeRaster(newDat,filename=files$combined[i],format='GTiff')
# #     rm(newDat)
# #     gc()
# #     print(paste0('Created file ',i,' of ',nrow(files)))
# #   } else {
# #     print(paste0('File ',i,' of ',nrow(files),' already exists'))
# #   }
# # } #Takes about 6 mins for 30 images
# 
# {a <- Sys.time()
#   sDat <- lapply(files$combined,function(datPath){
#     cNames <- names(brick(files$aqua[1])) #Channel names
#     sD <- brick(datPath) #Read data
#     names(sD) <- cNames #Add names
#     sD2 <- rasterToPoints(sD) #Convert to point dataframe
#     nMissing <- apply(sD2[,!grepl('(x|y)',colnames(sD2))],1,function(x) sum(is.na(x))) #Proportion of missing values in each row (cell)
#     # sD3 <- sD2[nMissing<5,] #Retains records with 4 or more data values
#     # if(nrow(sD3)==0) sD3 <- rbind(sD3,rep(NA,ncol(sD3)))
#     sD4 <- data.frame(doy=strsplit(x = datPath,split = c("(\\_|\\.)"))[[1]][2], #Day of year
#                       sD2,numNA=nMissing)
#     return(sD4)})
# Sys.time()-a} #Takes about 2 mins
# names(sDat) <- sapply(sDat,function(x) x$doy[1])
# 
# #Preview day 2
# sDat[[2]] %>% #sst on day 2
#   st_as_sf(coords=c('x','y')) %>% st_set_crs(4326) %>%
#   ggplot()+
#   geom_sf(aes(col=sst),size=0.1)+
#   geom_sf(data=dataBoundary,fill=NA)
# 
# sDat2 <- do.call('rbind',sDat) #Combine into single DF
# 
# # lwrLimits2 <- apply(dplyr::select(sDat2,chlor_a:sst),2,function(x) min(x,na.rm=TRUE)) #Minimum values from data
# # names(lwrLimits) <- c('chlor_a','nflh','poc','Rrs_412','Rrs_443','Rrs_469','Rrs_488','Rrs_531','Rrs_547','Rrs_555','Rrs_645','Rrs_667','Rrs_678','sst')
# 
# sDat2 %>% filter(numNA==0) %>% dplyr::select(chlor_a:sst) %>%
#   mutate(across(names(lwrLimits[lwrLimits<0]),~.x-(min(.x,na.rm=TRUE)*1.1))) %>%
#   slice_sample(n=10000) %>%
#   pivot_longer(everything()) %>%
#   mutate(value=log(value)) %>%
#   ggplot()+geom_density(aes(x=value))+facet_wrap(~name,scales='free')
# 
# sDat2 <- sDat2 %>%
#   # # #Rescales negative values be above 0.95*minimum positive value
#   # mutate(across(chlor_a:sst,~ifelse(.x<lwrLimits[names(lwrLimits)==cur_column()],lwrLimits[names(lwrLimits)==cur_column()]*0.95,.x)))
#   # Adds lowest value + 10% if lwrLimits < 0
#   mutate(across(names(lwrLimits[lwrLimits<0]),~.x-(min(.x,na.rm=TRUE)*1.1))) %>%
#   mutate(across(chlor_a:sst,log)) #Log-scale variables
# 
# sDatMat <- sDat2 %>% dplyr::select(chlor_a:sst) %>% as.matrix() #Data in matrix form
# sDat2$numNA <- apply(sDatMat,1,function(x) sum(is.na(x))) #Get number of NA values
# round(table(sDat2$numNA)/nrow(sDat2),4) #~50% cells missing all data, ~40% complete
# datMissing <- sDat2$numNA==13 #Locations missing all data
# 
# sDat2SomeMiss <- sDat2 %>% filter(numNA!=13,numNA!=0) #Some non missing
# sDat2NoMiss <- sDat2 %>% filter(numNA==0) #No missing
# 
# sDatMat <- sDat2SomeMiss %>% dplyr::select(chlor_a:sst) %>% as.matrix() #Partly missing data in matrix form
# 
# # n_components <- estim_ncpPCA(sDatMat, verbose = TRUE,method="Regularized",method.cv="gcv",ncp.min=1,ncp.max=13) #Takes about 5 mins
# # plot(1:13,n_components$criterion,xlab='Number of Dimensions',ylab='GCV Criterion',type='b',pch=19) #Looks like about 7 components is OK for prediction
# 
# #Same as imputePCA, but multicore. Splits X into slices of nPer each, passes to cluster cl, reassembles
# imputePCA_multi <- function(X,ncp=2,scale=TRUE,method='Regularized',nPer=10000,ncore=10){
#   # X <- sDatMat #Debugging
#   # nPer <- 10000
#   # ncore <- 10
#   # ncp <- 3
#   # scale <- TRUE
#   # method <- 'Regularized'
# 
#   #Assign rows to sample strata
#   sampNum <- c(rep(1:(nrow(X) %/% nPer),each=nPer), #Regular samples
#     rep((nrow(X) %/% nPer + 1),(nrow(X) %% nPer))) #"Leftovers"
# 
#   sampOrd <- sample(1:length(sampNum))
#   sampNum <- sampNum[sampOrd] #Randomize
# 
#   rNames <- rownames(X) #Save rownames
#   rownames(X) <- 1:nrow(X)
# 
#   IP <- function(x,ncp,scale,method) imputePCA(x,ncp,scale,method)$completeObs #Convenience function to pull only 1 matrix from output
# 
#   print('Breaking into separate matrices')
#   X <- lapply(1:max(sampNum),function(x) X[sampNum==x,]) #Break into separate matrices
#   if(ncore>1){
#     library(parallel)
#     cl <- makeCluster(ncore)
#     print(paste('Running imputePCA across',ncore,'clusters'))
#     clusterEvalQ(cl,library(missMDA)) #Loads missMDA on each cluster
#     for(i in 1:(ceiling(length(X)/ncore))){ #Break into sets. Tends to hang on large
#       use <- ((i-1)*ncore+1):pmin(length(X),(i*ncore)) #Which sets of matrices to use
#       X[use] <- parLapply(cl=cl,X=X[use],fun=IP,ncp=ncp,scale=scale,method=method) #Run imputePCA on separate matrices
#       cat('.') #Something to look at
#       # X <- parLapply(cl=cl,X=X,fun=IP,ncp=ncp,scale=scale,method=method) #Run imputePCA on separate matrices
#     }
#     stopCluster(cl); gc()
#   } else {
#     # X <- lapply(X,IP,ncp=ncp,scale=scale,method=method) #Run imputePCA on separate matrices
#     times <- rep(NA,length(X))
#     for(i in 1:length(X)){
#       a <- Sys.time()
#       X[[i]] <- IP(X[[i]],ncp=ncp,scale=scale,method=method)
#       b <- Sys.time()
#       times[i] <- difftime(b,a,units='secs')
#       print(paste0('Finished ',i,'. Time: ',round(times[i],2),'. Average time: ',round(mean(times,na.rm=TRUE),2),
#                    '. Estimated remaining (mins): ',round(mean(times,na.rm=TRUE),2)*(length(X)-i)/60 ))
#     }
#   }
#   X <- do.call(rbind,X) #Recombine into single matrix
# 
#   if(length(rNames)!=nrow(na.omit(X))) stop('NA values still present in dataframe')
# 
#   X <- X[order(as.numeric(rownames(X))),] #Reorder
#   rownames(X) <- rNames #Reset rownames
# 
#   print('Done')
#   return(X)
# }
# 
# #~3 mins
# sDatMat_imputed <- imputePCA_multi(sDatMat,ncp=5,scale=TRUE,method='Regularized',nPer = 10000,ncore = 15) #Impute missing data using 5 dimensions
# 
# sDat2SomeMiss[which(names(sDat2SomeMiss)=='chlor_a'):which(names(sDat2SomeMiss)=='sst')] <- sDatMat_imputed #Rejoin
# 
# rm(sDatMat_imputed); gc() #Cleanup
# 
# #Put back into a single dataframe
# sDat2 <- bind_rows(sDat2SomeMiss,sDat2NoMiss) %>% dplyr::select(-numNA) %>% arrange(doy,x,y)
# rm(sDat2SomeMiss,sDat2NoMiss,sDatMat,datMissing); gc()
# 
# sDatMat <- sDat2 %>% dplyr::select(chlor_a:sst) %>% as.matrix() #Data in matrix form again
# useThese <- unname(which(apply(sDatMat,1,function(x) !any(is.na(x))))) #Rows without missing values
# 
# #Decompose to PCs
# pca1 <- prcomp(sDatMat,scale. = TRUE,retx = FALSE)
# predict(pca1,sDatMat) %>% cor #Transformed values
# summary(pca1)
# 
# #Matrix to store PC values
# pcaMat <- matrix(NA,nrow = nrow(sDatMat), ncol=5,dimnames = list(rownames(sDatMat),paste0('PC',1:5)))
# pcaMat[useThese,] <- predict(pca1,sDatMat)[,1:ncol(pcaMat)] #Store PC values in matrix
# 
# sDat2 <- cbind(sDat2,pcaMat) #Combine PCs with sDat2
# 
# rm(sDatMat,pcaMat); gc()
# 
# #Add coordinate system - takes a minute or so
# sDat2 <- sDat2 %>% st_as_sf(coords=c('x','y')) %>%
#   dplyr::select(-chlor_a:-sst) %>% #Remove raw data, keep PCs
#   st_set_crs(4326) %>%
#   geom2cols(E,N,removeGeom=FALSE,epsg=3401) #AB UTM
# 
# Emean <- mean(unique(sDat2$E)); Nmean <- mean(unique(sDat2$N)) #Center values of E/N, used for scaling
# 
# sDat2 <- sDat2 %>%
#   mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
#   mutate(across(E:N,~round(.x))) %>%
#   st_transform(4326) %>% unite(loc,E,N) %>% #Back to WGS84
#   mutate(loc=as.numeric(factor(loc)))
# 
# withinBuff <- sDat2 %>% st_intersects(.,dataBoundary) %>% sapply(.,function(x) length(x)>0) #Points in sDat2 that are outside of the buffer
# sDat2 <- sDat2 %>% filter(withinBuff) #Filter out points outside of buffer
# 
# locs <- dplyr::select(sDat2,loc,sE,sN) %>% unique() #Unique locations in spectral data
# 
# # st_distance(locs[1:500,]) %>% data.frame() %>% dplyr::select(X1) %>%
# #   mutate(X1=as.numeric(X1)) %>%
# #   filter(X1!=0) %>% filter(X1==min(X1)) #Minimum nonzero distance is ~4km
# 
# rm(withinBuff,useThese); gc()
# 
# save(pca1,file='./data/PCAvals.Rdata') #pca values
# save(sDat,file = "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/sDat_raw.Rdata") #raw rasters
# #sf dataframe, plus other info
# save(sDat2,locs,Emean,Nmean,
#      file = "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/sDat2.Rdata")

load('./data/PCAvals.Rdata')
load("/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/sDat2.Rdata")

bottomWDat <- bottomWDat %>% mutate( #Match closest unique location in spectral data
  loc=apply(st_distance(bottomWDat,locs),1,which.min),
  minDist=apply(st_distance(bottomWDat,locs),1,min)) %>% 
  mutate(loc=ifelse(minDist>4000,NA,loc)) #Set nearest to NA if distance < 4 km (no matching cell)

# # PC1-3 changes through the year
# # NOTE: Many locations/times are completely missing
# sDat2 %>% 
#   filter(doy=='071'|doy=='141'|doy=='211') %>% 
#   pivot_longer(cols=c(PC1:PC5)) %>% 
#   ggplot()+
#   geom_sf(data=dataBoundary,fill=NA)+
#   geom_sf(aes(col=value),size=0.1)+
#   facet_grid(doy~name)

# ggplot()+geom_sf(dat=dataBoundary)+geom_sf(dat=locLookup)

# Fit PC models -----------------------------------

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
d=c(40,80,160,225), #Distances
N=c(30,25,20,8) #Points within each distance buffer
) %>%
  do.call('c',.)

#Knot locations
ggplot() + #Looks OK, but some bottom water points are outside boundary of data
  geom_sf(data=dataBoundary,fill=NA) +
  geom_sf(data=slice_sample(sDat2,n=5000),alpha=0.3)+
  geom_sf(data=bottomWDat,col='blue')+
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

#Breaks for 5 chunks of data
(doyBreaks <- sDat2 %>% slice(c(1,(nrow(sDat2) %/% 5)*c(1:5))) %>% pull(doy) %>% as.numeric)

sDat2 %>% st_drop_geometry() %>% group_by(doy) %>% 
  summarize(#propNA=sum(is.na(PC1))/n(),
            n=n()) %>% 
  mutate(doy=as.numeric(doy)) %>% 
  ggplot(aes(x=doy,y=n))+geom_point()+
  geom_vline(xintercept = doyBreaks)+labs(x='DOY',y='Number of points')

sDat2 <- sDat2 %>% #Add split data, convert DOY to numeric
  mutate(doySplit=cut(as.numeric(doy),breaks=doyBreaks,include.lowest=TRUE)) %>% 
  mutate(doy=as.numeric(doy))

sDatList <- lapply(levels(sDat2$doySplit),function(x){ #Divide into separate data chunks
  sDat2 %>% 
    st_drop_geometry() %>%
    mutate(across(c(sE,sN),unname)) %>% 
    filter(doySplit==x) %>% 
    dplyr::select(-doySplit)
})

#Per model: ~5 mins for 83 knots, 40 boundary loops, 10 layers for tensor product 
modForm <- "999~te(sN,sE,doy,bs=c('sf','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(40,10),d=c(2,1))+
  te(sN,sE,doy,bs=c('sw','tp'),xt=list(list(bnd=bound,nmax=100),NULL),k=c(40,10),d=c(2,1))"
#First te("sf") is for boundary, second te("sw") is for actual soap film

storageDir <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/models/"
# for(chunk in 1:(length(doyBreaks)-1)){ #For each chunk of data
#   for(pc in 1:5){ #For each PC
#     tempDat <- sDatList[[chunk]]
#     library(parallel)
#     cl <- makeCluster(15)
#     PCmod1 <- bam(formula(gsub('999',paste0('PC',pc),modForm)), knots=kts, data=tempDat, cluster=cl)  
#     stopCluster(cl)
#     
#     #Save coefficients
#     modList <- list(coefs=coef(PCmod1), #Coefficients
#                     vcv=vcov(PCmod1), #Covariance matrix
#                     smooths=PCmod1$smooth) #Smooth specifications
#     save(modList,file = paste0(storageDir,"modList",chunk,"_pc",pc,".Rdata")) 
#     
#     #Diagnostics
#     #https://stackoverflow.com/questions/45918662/plot-gam-from-mgcv-with-all-terms-true-within-a-function
#     png(paste0(storageDir,'summary',chunk,"_pc",pc,".png"),width=12,height=6,units='in',res=200)
#     # plot(PCmod1,scheme=2,too.far=0.01,pages=1,all.terms=TRUE)
#     plot(PCmod1,scheme=2,too.far=0.01,pages=1,all.terms=FALSE)
#     dev.off()
#     
#     #Residual plots
#     actual <- tempDat[,paste0('PC',pc)]
#     res <- residuals(PCmod1,type = 'pearson')
#     devRes <- residuals(PCmod1,type='deviance')
#     pred <- predict(PCmod1)
#     fit <- fitted(PCmod1)
#     if(length(devRes)>10000){
#       samp <- sample(1:length(devRes),10000) #Sub-samples to 10000
#       devRes <- devRes[samp]; pred <- pred[samp]; fit <- fit[samp]; res <- res[samp]; actual <- actual[samp]
#     } 
#     
#     #Residuals
#     png(paste0(storageDir,'residuals',chunk,"_pc",pc,".png"),width=8,height=8,units='in',res=200)
#     par(mfrow=c(2,2))
#     qqnorm(devRes,main='Deviance residuals Q-Q'); qqline(devRes)
#     plot(fit,res,xlab='Fitted',ylab='Resid',main='Fitted vs Residual'); abline(h=0,col='red')
#     hist(devRes,main='Deviance residuals')
#     plot(actual,pred,xlab='Actual',ylab='Predicted',main='Actual vs Predicted'); abline(0,1,col='red')
#     dev.off()
#     par(mfrow=c(1,1))
#     
#     capture.output({ #Model summary
#       print("SUMMARY---------------------------------")
#       summary(PCmod1) 
#     },file=paste0(storageDir,'results',chunk,"_pc",pc,".txt"))
#     # capture.output({ #GAM check results
#     #   print(" ")
#     #   print("GAM.CHECK-------------------------------")
#     #   k.check(PCmod1)# Runs into problems here. I don't think this work for soap film smoothers
#     # },file=paste0(storageDir,'results',chunk,"_pc",pc,".txt"),append=TRUE) 
#   }
# }


#OK, but not perfect. Looks like the spatial scale of the process is pretty small
#Turbulence patterns would probably have to have some kind of dynamic model

load(paste0(storageDir,"modList",'1',"_pc",'1',".Rdata"))

#PC1, predicted PC1, and residuals
sDatList[[1]] %>% dplyr::select(-sE,-sN) %>% 
  left_join(locs,by='loc') %>% 
  mutate(pred=predModList(modList,newdat = dplyr::select(sDatList[[1]],doy,sE,sN))) %>%
  mutate(resid=PC1-pred) %>% 
  filter(doy==66|doy==81|doy==99) %>%
  dplyr::select(PC1,pred,resid,doy,geometry) %>%
  pivot_longer(c(pred,resid)) %>%
  ggplot()+geom_sf(aes(geometry=geometry,col=value),size=0.1)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
  geom_sf(data=knotLocs)+geom_sf(data=dataBoundary,fill=NA)+
  facet_grid(name~doy)

#Residuals only
sDatList[[1]] %>% dplyr::select(-sE,-sN) %>% 
  left_join(locs,by='loc') %>% 
  mutate(pred=predModList(modList,newdat = dplyr::select(sDatList[[1]],doy,sE,sN))) %>%
  mutate(resid=PC1-pred) %>% 
  filter(doy==66|doy==81|doy==99) %>%
  ggplot()+geom_sf(aes(geometry=geometry,col=resid),size=0.1)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
  geom_sf(data=dataBoundary,fill=NA)+
  facet_wrap(~doy,ncol=1)

#mean(abs(residuals)): good idea for future maps

for(chunk in 1:length(sDatList)){
  for(pc in 1:5){
    load(paste0(storageDir,"modList",chunk,"_pc",pc,".Rdata"))
    p1 <- sDatList[[chunk]] %>% dplyr::select(-sE,-sN) %>% 
      mutate(pred=predModList(modList,newdat = dplyr::select(sDatList[[chunk]],doy,sE,sN))) %>% 
      mutate(resid=.data[[paste0('PC',pc)]]-pred) %>% 
      group_by(loc) %>% summarize(mae=mean(abs(resid))) %>% 
      ungroup() %>% left_join(locs,by='loc') %>% st_sf() %>% 
      ggplot()+
      geom_sf(aes(geometry=geometry,col=mae),size=0.6,shape=15)+
      geom_sf(data=knotLocs,col='red',size=0.5)+
      geom_sf(data=bottomWDat,size=0.5)+
      geom_sf(data=dataBoundary,fill=NA)+
      scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
      labs(title=paste0('errMap',chunk,"_pc",pc))
    
    ggsave(paste0(storageDir,'errMap',chunk,"_pc",pc,".png"),p1,scale = 2)
    rm(p1,modList); gc()
  }
}



#Get model performance stats
sResList <- lapply(sDatList,function(x){
  x$PC1 <- NA; x$PC2 <- NA; x$PC3 <- NA;
  x$PC4 <- NA; x$PC5 <- NA 
  return(x)
})

#Takes ~1 min
for(pc in 1:5){ 
  for(chunk in 1:(length(doyBreaks)-1)){
    load(paste0(storageDir,"modList",chunk,"_pc",pc,".Rdata"))
    nd <- sResList[[chunk]][,c('doy','sE','sN')]
    sResList[[chunk]][,paste0('PC',pc)] <- sDatList[[chunk]][,paste0('PC',pc)] - predModList(modList,newdat = nd) #Actual - predicted
    rm(modList,nd); gc() #Memory usage fairly heavy
  }
}

#Basic statistics
sResList %>% bind_rows(.id='chunk') %>% dplyr::select(contains('PC')) %>% 
  pivot_longer(everything()) %>% group_by(name) %>% summarize(mean=mean(value),mae=mae(value),rmse=rmse(value))

#Residuals over time
sResList %>% bind_rows(.id='chunk') %>% 
  group_by(doy) %>% 
  summarize(across(contains('PC'),list(med=median,upr=~quantile(.x,0.95),lwr=~quantile(.x,0.05)))) %>% 
  pivot_longer(contains('PC')) %>% 
  separate(name,c('name','stat'),sep='_') %>% 
  pivot_wider(names_from=stat,values_from = value) %>% 
  ggplot(aes(x=doy))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=med))+
  facet_wrap(~name,ncol=1)+
  geom_vline(xintercept = doyBreaks,linetype='dashed')+
  labs(x='Day of Year',y='Median Error')

#Residuals over space
sResList %>% bind_rows(.id='chunk') %>% 
  group_by(loc) %>% 
  summarize(across(contains('PC'),list(mae=mae,rmse=rmse))) %>% 
  pivot_longer(contains('PC')) %>% 
  separate(name,c('name','stat'),sep='_') %>% 
  left_join(locs,by='loc') %>% st_sf() %>% 
  ggplot()+
  geom_sf(aes(col=value),size=0.1)+
  # geom_sf(data=knotLocs)+ #Knot locations
  geom_sf(data=bottomWDat,alpha=0.3) + #Data locations
  geom_sf(data=dataBoundary,fill=NA)+
  facet_grid(name~stat)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)

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
l <- sort(unique(locs$loc)) #All unique locations
daylag <- 16 #Use 16-day lag

newDF <- expand_grid(doy=d,loc=l) %>% #Days/locations to predict at
  mutate(doy2=doy,doy=doy-daylag) %>% #doy = lagged day to use for predictions, doy 2 = Actual day of DO measurement
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN,-doySplit) %>% 
  left_join(st_drop_geometry(locs),by='loc') %>% 
  mutate(doySplit=cut(as.numeric(doy),breaks=doyBreaks,include.lowest=TRUE)) %>% 
  mutate(doySplit=as.numeric(doySplit))

#Separates into 2 df, one with values and one with NAs (quicker predictions)
newDF1 <- newDF %>% filter(!is.na(PC1)) %>% mutate(imputed=FALSE)  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% mutate(imputed=TRUE) #DF with NAs

for(pc in 1:5){ #Takes about 4 mins
  for(chunk in unique(newDF2$doySplit)){
    load(paste0(storageDir,"modList",chunk,"_pc",pc,".Rdata"))
    useThese <- newDF2$doySplit==chunk #indices to use
    nd <- newDF2 %>% filter(useThese) %>% dplyr::select(doy,sE,sN)
    newDF2[useThese,paste0('PC',pc)] <- predModList(modList,newdat = nd)
    rm(modList); gc() #Memory usage fairly heavy
  }
}

#Join back together
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
  left_join(locs,.,by='loc') %>% dplyr::select(-loc,-sE,-sN)

#Save a copy for YingJie
LLmapDat <- mapDat %>% 
  relocate(slice,minDO,meanDO,minDO_factor,meanDO_factor,propImputed) %>% 
  geom2cols()
save(LLmapDat,file = './data/LLmapDat.Rdata')

#Full copy for YingJie
LLmapDatFull <- newDF %>% left_join(.,locs,by='loc') %>% 
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
(p <- ggarrange(p1,p2,ncol=1))

ggsave(p,filename = './figures/mapLLpred.png',width=12,height=15)

# Get predictions from FR model ------------------------------

# load('./data/funRegMod.Rdata',verbose=TRUE) #Original FDA model uses 80 day range, but if we want to use 30 day lags:
load(file='./data/FRdat.Rdata')

fdat <- lapply(fdat,function(x){
  if(is.matrix(x)){
    return(x[,1:31])
  } else {
    return(x)
  }
})

basisType <- 'cr' #Cubic regression splines
#Fit FDA model 
bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+
                 s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+
                 s(dayMat,by=pcaMat4,bs=basisType)+
                 s(dayMat,by=pcaMat5,bs=basisType), 
               data=fdat)
summary(bWatMod)

dateRanges <- read.csv('./data/mappingDateRanges.csv') %>% #Get dates to predict on
  mutate(across(c(startDay,endDay), ~as.numeric(format(as.Date(paste0(.x,' 2014'),format='%B %d %Y'),format='%j')),.names='{.col}_doy'))
dateLabs <- with(dateRanges,paste0(startDay,' - ',endDay)) #Labels for date ranges
dateRanges <- with(dateRanges,data.frame(slice=rep(slice,nDays),doy=min(startDay_doy):max(endDay_doy)))

d <- dateRanges$doy #All days
l <- sort(unique(locs$loc)) #All unique locations
daylag <- 30 #30-day lag

#Dataframe of values to choose from
newDF <- expand_grid(doy=unique(do.call('c',lapply(d,function(x) x-(daylag:0)))),loc=l) %>% 
  left_join(st_drop_geometry(sDat2),by=c('doy','loc')) %>% dplyr::select(-sE,-sN) %>% 
  left_join(st_drop_geometry(locs),by='loc') %>% 
  mutate(doySplit=cut(as.numeric(doy),breaks=doyBreaks,include.lowest=TRUE)) %>% 
  mutate(doySplit=as.numeric(doySplit))

newDF1 <- newDF %>% filter(!is.na(PC1)) %>% mutate(imputed=FALSE)  #DF with values
newDF2 <- newDF %>% filter(is.na(PC1)) %>% mutate(imputed=TRUE) #DF with NAs

# cl <- makeCluster(15)
# {a <- Sys.time() #Takes about 4.6 mins for all 2.4 mil points
# newDF2 <- parLapply(cl=cl,X=list(PCmod1,PCmod2,PCmod3,PCmod4,PCmod5,PCmod6),fun=function(x,N){require(mgcv); predict.gam(x,newdata=N)},N=newDF2) %>% 
#   set_names(paste0('PC',1:6)) %>% bind_cols() %>% 
#   bind_cols(newDF2,.) %>% relocate(doy,loc,PC1:PC6,sE,sN) %>% mutate(imputed=TRUE)
# Sys.time()-a}; beep()
# stopCluster(cl)
# newDF <- bind_rows(newDF1,newDF2) %>% arrange(loc,doy)
# rm(newDF1,newDF2) #Cleanup


{ #Takes about 5.5 mins
  pb <- txtProgressBar(style=3)
  a <- Sys.time()
  pbi <- 1; pbl <- 5*length(unique(newDF2$doySplit))
  for(pc in 1:5){ 
    for(chunk in unique(newDF2$doySplit)){
      load(paste0(storageDir,"modList",chunk,"_pc",pc,".Rdata"))
      useThese <- newDF2$doySplit==chunk #indices to use
      nd <- newDF2 %>% filter(useThese) %>% dplyr::select(doy,sE,sN)
      newDF2[useThese,paste0('PC',pc)] <- predModList(modList,newdat = nd)
      rm(modList); gc() #Memory usage fairly heavy
      setTxtProgressBar(pb,pbi/pbl)
      pbi <- pbi+1
    }
  }
  b <- Sys.time()
  difftime(b,a)
}

#Join together
newDF <- bind_rows(newDF1,newDF2) %>% arrange(loc,doy)

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
}) %>% set_names(paste0('pcaMat',1:5))

datList <- c(datList,pcaMats); rm(pcaMats); gc()
datList$isImputed <- matrix(nd[,which(grepl('imputed',names(nd)))],ncol=ncol(datList$dayMat),byrow=TRUE) #Which values are imputed?

{a <- Sys.time() #Takes 4.6 mins to generate predictions for each day/location
  mapDat <- with(datList,data.frame(loc=loc,doy=doy,predDO=predict(bWatMod,newdata=datList))) 
  Sys.time()-a}

#Gets proportion of imputed data used in each prediction
mapDat$propImputed <- apply(datList$isImputed,1,mean)

#Full copy for YingJie
FRmapDatFull <- mapDat %>% left_join(.,locs,by='loc') %>% 
  dplyr::select(doy,propImputed,predDO,geometry) %>% st_as_sf() %>% 
  geom2cols()
save(FRmapDatFull,file = './data/FRmapDatFull.Rdata')

mapDat <- mapDat %>% 
  mutate(slice=as.numeric(cut(doy,breaks=with(dateRanges,c(startDay_doy[1],endDay_doy)),include.lowest=TRUE))) %>% 
  left_join(.,dateRanges,by='slice') %>% 
  group_by(slice,loc) %>% summarize(minDO=min(predDO),meanDO=mean(predDO),propImputed=mean(propImputed)) %>% 
  ungroup() %>% 
  mutate(across(c(minDO,meanDO),~cut(.x,breaks=c(min(.x,na.rm=TRUE),1,2,3,4,5,max(.x,na.rm=TRUE)),labels=c('<1','1-2','2-3','3-4','4-5','>5'),include.lowest=TRUE),
                .names='{.col}_factor')) %>% 
  mutate(slice=factor(slice,lab=dateLabs)) %>% 
  left_join(locs,.,by='loc') %>% dplyr::select(-loc,-sE,-sN)

#Save data for Yingie
FRmapDat <- mapDat %>% 
  relocate(slice,minDO,meanDO,minDO_factor,meanDO_factor,propImputed) %>% 
  mutate(propImputed=NA) %>% 
  geom2cols()
save(FRmapDat,file = './data/FRmapDat.Rdata')

#Make figure
(p1 <- ggplot(mapDat)+
  geom_sf(aes(col=minDO_factor,geometry=geometry),size=0.2)+
  facet_wrap(~slice,ncol=3)+
  scale_colour_brewer(type='seq',palette = 'RdBu')+
  labs(title='Functional Regression - min(DO)',col='DO (mg/L)'))
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


