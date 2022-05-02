#THIRD APPROACH TO ANALYSIS: MODEL BOTTOM DO BASED ON RAW SATELLITE CHANNELS + GAP-FILLING
#WRITTEN BY SR, SPRING 2021

library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(ggmap)
library(ggpubr)
library(animation)
library(mgcv)
library(vegan)
library(gstat)

setwd("~/Documents/hypoxiaMapping")

source('helperFunctions.R')

#Load data from saved file
# load('./data/all2014_2.Rdata') #November 2021

#GAM imputation ------------------------------------------

#Load smoothers
# load('./data/PCmods.RData') #Thin-plate splines
# load('./data/PCmodsSoap.RData') #Soap film - this uses the smoothers fit to the data at each sampling location, not all the spectral data
load('./data/all2014_3.Rdata') #Newer - March 2022, sDat here contains only spectral data from water-testing locations
load('./data/PCAvals.Rdata')
load("/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/sDat2.Rdata") #All spectral data

# > Emean #Currently
# [1] 2277503
# > Nmean
# [1] 3382200
# > Emean #Should be...
# [1] 2570555
# > Nmean
# [1] 3329257

storageDir <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/models/"
doyBreaks <- sDat2 %>% slice(c(1,(nrow(sDat2) %/% 5)*c(1:5))) %>% pull(doy) %>% as.numeric

bottomWDat <- bottomWDat %>% mutate( #Match closest unique location in spectral data
  loc=apply(st_distance(bottomWDat,locs),1,which.min),
  minDist=apply(st_distance(bottomWDat,locs),1,min)) %>% 
  mutate(loc=ifelse(minDist>4000,NA,loc)) %>% #Set nearest to NA if distance < 4 km (no matching cell)
  filter(!is.na(loc)) %>% 
  mutate(chunk=as.numeric(cut(as.numeric(doy),breaks=doyBreaks,label=1:5,include.lowest=TRUE))) %>% 
  select(-E,-N,-sE,-sN) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>%
  mutate(sE=(E-Emean)/1000,sN=(N-Nmean)/1000) %>% #Center E/N and convert to km
  st_transform(4326) #Back to WGS84

sDat2 <- sDat2 %>% mutate(chunk=as.numeric(cut(as.numeric(doy),breaks=doyBreaks,label=1:5,include.lowest=TRUE)))

# #Check alignment
# sDat2 %>% slice_sample(n=10000) %>%
#   ggplot()+geom_sf()+
#   geom_sf(data=bottomWDat,col='red')
# 
# sDat2 %>% slice_sample(n=10000) %>%
#   geom2cols(E2,N2,removeGeom=FALSE,epsg=3401) %>%
#   ggplot(aes(x=E2,y=N2))+geom_point()+
#   geom_point(data=geom2cols(bottomWDat,E2,N2,removeGeom=FALSE,epsg=3401),col='red')
# 
# sDat2 %>% slice_sample(n=10000) %>%
#   ggplot(aes(x=sE,y=sN))+geom_point()+
#   geom_point(data=bottomWDat,col='red')

#Get predictions of PCs at all locations through the entire season. If PC values missing, fill using PC model
lags <- 0:80 #Use lag days of 0-80

#Storage matrices for PC data
sDatList <- lapply(1:5,function(x) matrix(NA,nrow=nrow(bottomWDat),ncol=length(lags),dimnames=list(bottomWDat$YEID,paste0('L',lags)))) %>% 
  set_names(nm=paste0('PC',1:5))

{
  pb <- txtProgressBar(style=3)
  for(i in 1:nrow(bottomWDat)){ #Gets existing data
    if(bottomWDat$loc[i] %in% sDat2$loc){
      sDat2Loc <- st_drop_geometry(sDat2[sDat2$loc==bottomWDat$loc[i],])
      sDat2Loc$lag <- bottomWDat$doy[i]-as.numeric(sDat2Loc$doy)
      sDat2Loc <- sDat2Loc[sDat2Loc$lag %in% lags,]
      if(nrow(sDat2Loc)!=0){
        for(pc in 1:5){ #For each PC
          for(lagDay in 1:nrow(sDat2Loc)){
            sDatList[[pc]][i,sDat2Loc$lag[lagDay]] <- sDat2Loc[lagDay,paste0('PC',pc)]
          }
        }
      }
    } 
    setTxtProgressBar(pb,i/nrow(bottomWDat))
  }
  close(pb)
}

imputed <- is.na(sDatList[[1]]) #Are variables imputed?

{
  pb <- txtProgressBar(style=3)
  for(ch in 1:(length(doyBreaks)-1)){ #For each chunk
    doyRange <- doyBreaks[c(ch,ch+1)] #Range of days used in this chunk
    modPaths <- dir(storageDir,pattern = paste0('modList',ch),full.names = TRUE) #Paths to models
    modList <- lapply(modPaths,function(x){load(x); return(modList)}) #Load models into list
    
    dayMat <- t(sapply(bottomWDat$doy,function(x) x-lags)) #Matrix of days matching lag day values for each location
    eMat <- outer(bottomWDat$sE,rep(1,length(lags))) #Matrix of E and N indices
    nMat <- outer(bottomWDat$sN,rep(1,length(lags)))
    useMat <- dayMat>=doyRange[1] & dayMat<doyRange[2] & is.na(sDatList[[1]]) #Matrix of days/locations to get values for, excluding existing values
    predDF <- data.frame(doy=dayMat[useMat],sE=eMat[useMat],sN=nMat[useMat]) #DF to use for prediction
    
    predDF <- lapply(modList,predModList,newdat=predDF) %>% #Get predictions
      set_names(nm = paste0('PC',1:5)) %>% 
      bind_cols() %>% bind_cols(predDF,.)
    for(pc in 1:5) sDatList[[pc]][useMat] <- predDF[,paste0('PC',pc)]
    setTxtProgressBar(pb,ch/(length(doyBreaks)-1))
  }
  close(pb)
}
sDatList[[1]][1:10,1:10] #Works

sDat <- sDatList[[1]] %>% as.data.frame() %>% rownames_to_column('YEID') %>% 
  mutate(doy=bottomWDat$doy) %>% 
  st_sf(geometry=bottomWDat$geometry) %>% 
  pivot_longer(c(-YEID,-doy,-geometry)) %>% rename(lag=name) %>% 
  mutate(lag=as.numeric(gsub('L','',lag)))

sDat <- imputed %>% as.data.frame() %>% rownames_to_column('YEID') %>% 
  pivot_longer(c(-YEID),values_to = 'imputed') %>% rename(lag=name) %>% 
  mutate(lag=as.numeric(gsub('L','',lag))) %>% select(imputed) %>% 
  bind_cols(sDat,.)

sDat <- lapply(sDatList,function(x){
  x %>% as.data.frame() %>% rownames_to_column('YEID') %>% 
    pivot_longer(-YEID) %>% rename(lag=name) %>% 
    mutate(lag=as.numeric(gsub('L','',lag))) %>% 
    pull(value)}) %>% bind_cols(sDat,.)

sDat <- sDat %>% mutate(doy=doy-lag) %>% select(-lag) %>% 
  mutate(date=as.Date(paste('2014',round(doy),sep='-'),format='%Y-%j'))
  
#Fit model of DO to gap-filled PCs ----------------------------------

lags <- 0:80 #Try 0 to 80 day lags

fitLagMods <- function(i,dat=sDat,interaction=FALSE){
  #Copy of spectral data
  sDatTemp <- st_drop_geometry(dat) %>% 
    select(YEID,doy,contains('PC')) %>% 
    mutate(doy=doy+i) %>% #Add to go back, subtract to go forward
    unite(ID,c('YEID','doy'),sep='-')
  
  bWatTemp <- bottomWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% 
    left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1))
  
  if(!interaction){
    #Fit simple linear models use PC1:5
    bMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=bWatTemp)
    
  } else {
    #Interaction models
    bMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=bWatTemp)
  }
  if(any(is.na(coef(bMod)))|length(coef(bMod))>=nrow(bWatTemp)-5){
    warning('NA coefs or more coefs than data. Model fit singular.')
    return(NA)
  }
  return(bMod)
}

modList1 <- lapply(lags,fitLagMods,interaction=FALSE)
modList2 <- lapply(lags,fitLagMods,interaction=TRUE)

#Get plots of MSE, MAE, and R-squared
(p <- list(lag=lags,
     bind_rows(lapply(modList1,function(x) data.frame(rmse_simple=rmse(x),mae_simple=mae(x),r2_simple=getR2(x)))),
     bind_rows(lapply(modList2,function(x) data.frame(rmse_interaction=rmse(x),mae_interaction=mae(x),r2_interaction=getR2(x))))) %>% 
  bind_cols() %>% 
  pivot_longer(-lag) %>% separate(name,c('stat','modType'),sep='_') %>% 
  mutate(stat=factor(stat,levels=c('mae','rmse','r2'),labels=c('Mean Absolute Error','Root Mean Square Error','R Squared'))) %>%
  mutate(modType=factor(modType,levels=c('simple','interaction'),labels = c('Simple','Interaction'))) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~stat,scales = 'free_y',ncol=1, strip.position = 'right')+
  labs(x='Time lag',y=NULL,col='Model Type'))
ggsave('./figures/lagPCAmod_gapfill.png',p,width=8,height=8)

#Global minimum at 80 day lag, 54 for interaction model. This seems unlikely, but using this for now
#Simple model
(bestDay <- lags[which.min(sapply(modList1,mae)[1:31])]) #Minimum mae occurs on day 80; if using only first 30 days, minimum on day 16
lags[which.min(sapply(modList1,rmse)[1:31])] #Minimum rmse occurs on day 80; 18 if using first 30 days

# #Interaction model
# (bestDay <- which.min(sapply(modList2,mae))) #Minimum mae occurs on day 53
# which.min(sapply(modList2,rmse)) #Day 53

m1 <- modList1[[bestDay+1]] #Save model from that day
save(m1,file = './data/lagLinMod.Rdata')

# Cross-validation (use 70%, predict on 30%)

#Function to fit models on random 70% of data, return prediction differences for remaining 30%
fitLagModsCV <- function(i,dat=sDat,interaction=FALSE,use=0.7){

  require(tidyverse)
  require(sf)
  
  #Copy of spectral data
  sDatTemp <- st_drop_geometry(dat) %>% 
    select(YEID,doy,contains('PC')) %>% 
    mutate(doy=doy+i) %>% #Add to go back, subtract to go forward
    unite(ID,c('YEID','doy'),sep='-')
  
  # sWatTemp <- surfWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% #Copies of water data
  #   left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1)) %>% mutate(predSet=makeTF(.,1-use))
    
  bWatTemp <- bottomWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% 
    left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1)) %>% mutate(predSet=makeTF(.,1-use))
  
  # sWatPredSet <- sWatTemp %>% filter(predSet) #Take only samples from prediction set
  bWatPredSet <- bWatTemp %>% filter(predSet) 
  
  # sWatTemp <- sWatTemp %>% filter(!predSet) #Remove prediction set
  bWatTemp <- bWatTemp %>% filter(!predSet)
  
  if(!interaction){
    #Fit simple linear models use PC1:5
    # sMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=sWatTemp)
    bMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=bWatTemp)
    
  } else {
    #Interaction models
    # sMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=sWatTemp)
    bMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=bWatTemp)
  }
  # if(any(is.na(coef(sMod)))|length(coef(sMod))>=nrow(sWatTemp)-5){
  #   warning('NA coefs or more coefs than data. Model fit singular.')
  #   return(list(surface=NA,bottom=NA))
  # }
  
  #Predict on withheld dataset, get difference
  # sModPreds <- sWatPredSet %>% mutate(pred=predict(sMod,newdata=sWatPredSet),diff=pred-DO)
  bModPreds <- bWatPredSet %>% mutate(pred=predict(bMod,newdata=bWatPredSet),diff=pred-DO)
  
  return(list(
    # surfaceDiff=sModPreds$diff,
    bottomDiff=bModPreds$diff,
    #R2 (not adjusted, but all models have the same # of coefs) - returning negatives for some reason
    # bottomR2=1-with(bModPreds,sum(diff^2)/sum((DO-mean(DO))^2)),
    bottomR2=summary(lm(DO~pred,data=bModPreds))$r.squared
    ))
}

#Function to run this in parallel, return as a dataframe. (runs 1000 replications)
cvPredErrs <- function(i,inter=FALSE){
  source('helperFunctions.R')
  lapply(replicate(n=1000,expr=fitLagModsCV(i,interaction = inter),simplify = FALSE),function(x){
    ret <- data.frame(rmseBottom=rmse(x$bottomDiff),maeBottom=mae(x$bottomDiff),r2Bottom=x$bottomR2)
    return(ret)
  })
}

# # Get prediction errors for surface/bottom, using RMSE, MAE, and R2
# library(parallel) #Takes about 10-15 mins
# cluster <- makeCluster(10)
# clusterExport(cluster,c('fitLagModsCV','sDat','bottomWDat'))
# cvPredList <- parLapply(cl=cluster,lags,cvPredErrs) #No interactions
# cvPredList2 <- parLapply(cl=cluster,lags,cvPredErrs,inter=TRUE) #Interactions
# stopCluster(cluster)
# beepr::beep(1)
# cvPredList <- lapply(cvPredList,function(x) do.call('rbind',x))
# cvPredList2 <- lapply(cvPredList2,function(x) do.call('rbind',x))
# 
# cvPredList <- cvPredList %>% bind_rows(.id='lag') %>% mutate(lag=as.numeric(lag)) %>% pivot_longer(-lag) %>%
#   mutate(errType=factor(name,labels=c('MAE','R2','RMSE'))) %>% select(-name) %>% mutate(modType='Simple')
# 
# cvPredList2 <- cvPredList2 %>% bind_rows(.id='lag') %>% mutate(lag=as.numeric(lag)) %>% pivot_longer(-lag) %>%
#   mutate(errType=factor(name,labels=c('MAE','R2','RMSE'))) %>% select(-name) %>% mutate(modType='Interaction')
# save(cvPredList,cvPredList2,file='./data/cvPredLists.Rdata')
load('./data/cvPredLists.Rdata')

# bind_rows(cvPredList,cvPredList2) %>%  #Interaction model clearly has something wrong with it (probaby overfitting)
#   ggplot(aes(x=lag,y=value,col=modType))+
#   geom_point(alpha=0.1,position=position_dodge(width=0.5))+
#   facet_wrap(~errType)+
#   labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged linear Model',col='Model Type')+
#   coord_cartesian(ylim=c(0,5))

# bind_rows(cvPredList,cvPredList2) %>% 
#   group_by(lag,errType,modType) %>% 
#   summarize(med=mean(value)) %>% ungroup() %>% 
#   pivot_wider(names_from=errType,values_from = c(med)) %>% 
#   group_by(modType) %>% 
#   mutate(minMAE=MAE==min(MAE),minRMSE=RMSE==min(RMSE)) %>% 
#   filter(minMAE|minRMSE) %>% select(-contains('min')) %>% 
#   data.frame() %>% 
#   relocate(lag,modType,RMSE,MAE,R2) 

#Within-sample error
WSE <- bind_cols(lag=lags,
  bind_rows(lapply(modList1,function(x) data.frame(MAE=mae(x),R2=getR2(x),RMSE=rmse(x))))  
) %>% pivot_longer(-lag,names_to = 'errType', values_to = 'mean_wse')

#Out of sample error
OOSE <- cvPredList %>% 
  # group_by(lag,errType) %>% 
  # summarize(mean=mean(value),med=median(value),max=max(value),min=min(value),iqr=IQR(value)) %>% 
  # summarize(med_oose=median(value),max_oose=max(value),min_oose=min(value)) %>% 
  mutate(lag=lag-1) %>% ungroup()

#Arrange for Yingjie
WSE <- WSE %>% transmute(lag,errType,error=mean_wse,errType2='Within Sample')
OOSE <- OOSE %>% transmute(lag,errType,error=value,errType2='Out of Sample')

#Join together
llErrLag <- bind_rows(WSE,OOSE)
save(llErrLag,file='./data/llErrLag.Rdata')

llErrLag %>% group_by(lag,errType,errType2) %>% 
  summarize(medErr=median(error),maxErr=max(error),minErr=min(error)) %>% 
  ggplot(aes(x=lag))+
  geom_ribbon(aes(ymax=maxErr,ymin=minErr,fill=errType2),alpha=0.3)+
  geom_line(aes(y=medErr,col=errType2))+
  facet_wrap(~errType,ncol=1,scales='free_y')+
  labs(x='Time Lag',y='Error',col='Error Type',fill='Error Type')+
  scale_colour_manual(values=c('grey10','red'))+
  scale_fill_manual(values=c('grey10','red'))

# Fit FR model of DO to gap-filled PCs ------------------------------------

NdayLag <- 80 #80 days in past
NdayForward <- 0 #0 days in future
dayLags <- -NdayForward:NdayLag

#Data for functional regression
fdat <- list(DO_bottom=bottomWDat$DO,
             # DO_surf=surfWDat$DO,
             dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
             pcaMat1=sDatList[[1]],pcaMat2=sDatList[[2]],
             pcaMat3=sDatList[[3]],pcaMat4=sDatList[[4]],
             pcaMat5=sDatList[[5]],
             doy=bottomWDat$doy,
             sE=bottomWDat$sE,sN=bottomWDat$sN,
             maxDepth=bottomWDat$maxDepth)

sum(imputed)/prod(dim(imputed)) #77% imputed
plot(jitter(bottomWDat$doy),jitter(apply(imputed,1,mean)), #More imputation towards the end of the season
     xlab='DOY',ylab='Prop imputed',pch=19,cex=0.5)
save(fdat,file='./data/FRdat.Rdata')

basisType <- 'cr' #Cubic regression splines
#Fit FDA models 
bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+
                 s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+
                 s(dayMat,by=pcaMat4,bs=basisType)+
                 s(dayMat,by=pcaMat5,bs=basisType), 
               data=fdat) #Bottom water
save(bWatMod,file = './data/funRegMod.Rdata')

summary(bWatMod) #R-squared of about 0.63
par(mfrow=c(2,2)); gam.check(bWatMod); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad. Gets slightly worse at low fitted values
plot(bWatMod,scheme=1,pages=1) #Interesting cyclical patterns, esp for pc 4 - 5

#Use smoothPred to get FR plots from each smoother
pvals <- unname(round(summary(bWatMod)$s.table[,4],3))
pvals <- ifelse(pvals==0,'<0.001',paste0('=',as.character(pvals)))

(p1 <- lapply(1:5,function(i){
  d <- expand.grid(dayMat=0:80,p=1) #Dataframe
  names(d)[2] <- paste0('pcaMat',i) #Change name of by variable
  smoothPred(m=bWatMod,dat=d,whichSmooth=i)}) %>% 
  set_names(paste0('PCA',1:5,' (p',pvals,')')) %>% 
  bind_rows(.id='PC') %>% 
  select(-contains('pcaMat')) %>% 
  ggplot(aes(x=dayMat))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred))+facet_wrap(~PC)+geom_hline(yintercept=0,col='red',linetype='dashed')+
  labs(x='Day (lag)',y='Effect'))

(p2 <- data.frame(pred=predict(bWatMod),actual=fdat$DO_bottom) %>% 
  ggplot()+geom_point(aes(x=pred,y=actual))+
  geom_abline(intercept = 0, slope = 1)+
  labs(x='Predicted Bottom DO',y='Actual Bottom DO'))

(p <- ggarrange(p1,p2,ncol=2))
ggsave('./figures/frPCA_gapfill.png',p,width=10,height=5)
ggsave('./figures/frPCA_gapfill2.png',p1,width=10,height=5)

# Cross-validation (use 70%, predict on 30%)

#Function to fit models on random 70% of data, return prediction differences for remaining 30%
#uses fDat from above, and samples within it
fitFRmodCV <- function(i,dat=fdat,use=0.7){ #i doesn't do anything, just used for parallel processing

  require(mgcv)
  source('helperFunctions.R')
  
  #Copy of spectral data
  predSet <- makeTF(data.frame(dat$DO_bottom),(1-use))
  
  getRows <- function(x,choose) if(any(class(x)=='matrix')) return(x[choose,]) else return(x[choose])
  
  predDat <- lapply(fdat,getRows,choose=predSet) #Data to predict on
  dat <- lapply(fdat,getRows,choose=!predSet) #Data to fit with
  
  basisType <- 'cr' #Cubic regression splines
  #Fit FDA models 
  bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+
                   s(dayMat,by=pcaMat2,bs=basisType)+
                   s(dayMat,by=pcaMat3,bs=basisType)+
                   s(dayMat,by=pcaMat4,bs=basisType)+
                   s(dayMat,by=pcaMat5,bs=basisType),
                 data=dat) #Bottom water
  
  #Predict on withheld dataset, get difference
  predDat$pred <- predict(bWatMod,newdata = predDat)
  predDat$diff <- predDat$pred-predDat$DO_bottom
  
  ret <- with(predDat,c(rmse(diff),mae(diff),summary(lm(DO_bottom~pred,data=predDat))$r.squared))
  names(ret) <- c('RMSE','MAE','R2')
  return(ret)
}

library(parallel)
cluster <- makeCluster(10)
clusterExport(cluster,c('fitFRmodCV','fdat'))
cvPredList3 <- parLapply(cl=cluster,1:1000,fitFRmodCV)  #Takes only a few seconds to run
stopCluster(cluster)

pivot_longer(bind_rows(cvPredList3),RMSE:R2) %>% group_by(name) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value),iqr=IQR(value))

#What is the shortest lag time that we could use for FR and get similar results?

#Fit FDA models with different lags (10 days - 80 days)
lags <- 10:80

#Function to get error from model of lagged data. default = within-sample, 0<test<1 = out of sample
getLagErr <- function(l,f,test=0,N=30){
  library(mgcv)
  f$dayMat <- f$dayMat[,1:l]
  f$pcaMat1 <- f$pcaMat1[,1:l]
  f$pcaMat2 <- f$pcaMat2[,1:l]
  f$pcaMat3 <- f$pcaMat3[,1:l]
  f$pcaMat4 <- f$pcaMat4[,1:l]
  f$pcaMat5 <- f$pcaMat5[,1:l]
  
  basisType <- 'cr' #Cubic regression splines
  
  if(test==0){
    #Fit FDA models
    b <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+
               s(dayMat,by=pcaMat2,bs=basisType)+
               s(dayMat,by=pcaMat3,bs=basisType)+
               s(dayMat,by=pcaMat4,bs=basisType)+
               s(dayMat,by=pcaMat5,bs=basisType),
             data=f)
    
    ret <- c(rmse(b),mae(b),summary(b)$r.sq)
    names(ret) <- c('RMSE','MAE','R2')
    return(ret)
  } else {
    t(replicate(N,{
      #Get testing/training data indices
      ndat <- length(f$DO_bottom)
      testThese <- sort(sample(1:ndat,round(test*ndat)))
      trainThese <- (1:ndat)[!(1:ndat %in% sample(1:ndat,round(test*ndat)))]
      
      selectDat <- function(x,i){
        if(class(x)=='numeric') x <- x[i] else x <- x[i,]
      }
      
      testdat <- lapply(f,selectDat,i=testThese)
      traindat <- lapply(f,selectDat,i=trainThese)
      
      #Fit model
      b <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+
                 s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+
                 s(dayMat,by=pcaMat4,bs=basisType)+
                 s(dayMat,by=pcaMat5,bs=basisType),
               data=traindat)
      
      #Get difference between predicted/actual from test data
      res <- testdat$DO_bottom-predict(b,newdata=testdat)
      ret <- c(rmse(res),mae(res),summary(b)$r.sq)
      names(ret) <- c('RMSE','MAE','R2')
      return(ret)
    }))
  }
}

WSE <- lapply(lags,getLagErr,f=fdat) %>% #Training data (within-sample) error - takes a few seconds
  bind_rows() %>% mutate(lag=lags) %>%
  pivot_longer(RMSE:R2,names_to='errType',values_to='error') %>% 
  mutate(errType2='Within Sample')

# cluster <- makeCluster(10)
# clusterExport(cluster,c('getLagErr','rmse','mae'))
# a <- Sys.time()
# OOSE <- parLapply(cl=cluster,X = lags, fun = getLagErr,f=fdat,test=0.3,N=1000) #Testing data (out of sample) error - takes 17 mins
# stopCluster(cluster)
# b <- Sys.time()
# difftime(b,a)
# save(OOSE,file='./data/cvPredLists2.Rdata')
load('./data/cvPredLists2.Rdata')

OOSE <- OOSE %>% lapply(.,data.frame) %>% set_names(nm=as.character(lags)) %>%
  bind_rows(.id = 'lag') %>%
  mutate(lag=as.numeric(lag)) %>%
  pivot_longer(RMSE:R2,names_to='errType',values_to = 'error') %>% 
  mutate(errType2='Out of Sample') 

#Join together
fdaErrLag <- bind_rows(WSE,OOSE)
save(fdaErrLag,file='./data/fdaErrLag.Rdata')

fdaErrLag %>% group_by(lag,errType,errType2) %>% 
  summarize(medErr=median(error),maxErr=max(error),minErr=min(error)) %>% 
  ggplot(aes(x=lag))+
  geom_ribbon(aes(ymax=maxErr,ymin=minErr,fill=errType2),alpha=0.3)+
  geom_line(aes(y=medErr,col=errType2))+
  facet_wrap(~errType,ncol=1,scales='free_y')+
  labs(x='Time Lag',y='Error',col='Error Type',fill='Error Type')+
  scale_colour_manual(values=c('grey10','red'))+
  scale_fill_manual(values=c('grey10','red'))

#Conclusion: Data from 30 days previous is better, but not a huge amount better than, data from a 10-day span

#Compare within-sample performance of models --------------------------

#RMSE
min(sapply(modList1,rmse)) #Lagged linear regression - PCA
min(sapply(modList2,rmse)) #Lagged linear regression - PCA with interactions
rmse(bWatMod) #Functional regression - PCA

#Which days do these occur on?
lags[which.min(sapply(modList1,rmse))] #12-day lag
lags[which.min(sapply(modList2,rmse))] #0-day lag

#MAE
min(sapply(modList1,mae)) #Lagged linear regression - PCA
min(sapply(modList2,mae)) #Lagged linear regression - PCA with interactions
mae(bWatMod) #Functional regression - PCA

#R2
max(sapply(modList1,getR2)) #Lagged linear regression - PCA
max(sapply(modList2,getR2),na.rm=TRUE) #Lagged linear regression - PCA
summary(bWatMod)$r.sq #Functional regression - PCA

#df
sapply(modList1,getDF)[which.min(sapply(modList1,rmse))]
sapply(modList2,getDF)[which.min(sapply(modList2,rmse))]
bWatMod$df.residual

#Compare within-sample to out-of-sample error, using only simple lagged-linear model
load('./data/llErrLag.Rdata')
load('./data/fdaErrLag.Rdata')
errLag <- bind_rows(
  mutate(modType='Lagged Linear',llErrLag),
  mutate(modType='Functional Data Analysis',fdaErrLag)
)
NOTE <- 'Data are quite different from each other, so had to save as separate objects. llModErr = Out-of-sample and within-sample error for lagged linear models. fdaErr = Out-of-sample and within-sample error for 500 random draws of FDA model.'
save(errLag,NOTE,file = './data/errDat.Rdata')

load('./data/errDat.Rdata')
(p <- errLag %>% 
    group_by(lag,errType,errType2,modType) %>% 
    summarize(medErr=median(error),maxErr=max(error),minErr=min(error)) %>% 
    ggplot(aes(x=lag))+
    geom_ribbon(aes(ymax=maxErr,ymin=minErr,fill=errType2),alpha=0.3)+
    geom_line(aes(y=medErr,col=errType2))+
    facet_grid(errType~modType,scales='free_y')+
    labs(x='Time lag',y='Model Accuracy')+
    labs(x='Time Lag',y='Error',col='Error Type',fill='Error Type')+
    scale_colour_manual(values=c('grey10','red'))+
    scale_fill_manual(values=c('grey10','red'))
  )
ggsave('./figures/outOfSamp_comparison.png',p,width=8,height=8)

