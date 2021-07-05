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
library(missMDA)
library(gstat)

setwd("~/Documents/hypoxiaMapping")

source('helperFunctions.R')

#Load data from saved file
load('./data/all2014.Rdata')

#GAM imputation ------------------------------------------

#Load smoothers
load('./data/PCmods.RData')

#Get predictions of PCs at all locations through the entire season. If PC values missing, fill using PC model
sDat <- sDat %>% mutate(predPC1=predict(PCmod1,newdata=.),predPC2=predict(PCmod2,newdata=.)) %>% 
  mutate(predPC3=predict(PCmod3,newdata=.),predPC4=predict(PCmod4,newdata=.),predPC5=predict(PCmod5,newdata=.)) %>% 
  mutate(PC1=ifelse(gap,predPC1,PC1),PC2=ifelse(gap,predPC2,PC2),PC3=ifelse(gap,predPC3,PC3),PC4=ifelse(gap,predPC4,PC4),PC5=ifelse(gap,predPC4,PC5)) %>% 
  select(-predPC1:-predPC5) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j')) 

#Fit model of DO to gap-filled PCs ----------------------------------

lags <- 0:30 #Try 0 to 30 day lags
fitLagMods <- function(i,dat=sDat,interaction=FALSE){
  #Copy of spectral data
  sDatTemp <- st_drop_geometry(dat) %>% 
    select(YEID,doy,contains('PC')) %>% 
    mutate(doy=doy+i) %>% #Add to go back, subtract to go forward
    unite(ID,c('YEID','doy'),sep='-')
  
  sWatTemp <- surfWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% #Copies of water data
    left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1))
  bWatTemp <- bottomWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% 
    left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1))
  
  
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
p1 <- data.frame(lag=lags,
                 surface1=sapply(modList1,function(i) mae(i$surface)),
                 bottom1=sapply(modList1,function(i) mae(i$bottom)),
                 surface2=sapply(modList2,function(i) mae(i$surface)),
                 bottom2=sapply(modList2,function(i) mae(i$bottom))) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='Mean Absolute Error',col='Model Type')

#Get plots of MSE and RMSE
p2 <- data.frame(lag=lags,
                 surface1=sapply(modList1,function(i) rmse(i$surface)),
                 bottom1=sapply(modList1,function(i) rmse(i$bottom)),
                 surface2=sapply(modList2,function(i) rmse(i$surface)),
                 bottom2=sapply(modList2,function(i) rmse(i$bottom))) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='Root mean squared error',col='Model Type')

#Get plots of MSE and R2
p3 <- data.frame(lag=lags,
                 surface1=sapply(modList1,function(i) getR2(i$surface)),
                 bottom1=sapply(modList1,function(i) getR2(i$surface)),
                 surface2=sapply(modList2,function(i) getR2(i$surface)),
                 bottom2=sapply(modList2,function(i) getR2(i$surface))) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='R-squared',col='Model Type')

(p <- ggarrange(p1,p2,p3,ncol=1,common.legend=TRUE,legend='bottom') )
ggsave('./figures/lagPCAmod_gapfill.png',p,width=8,height=8)

#Similar to models using raw data: no local minimum appears

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
  
  sWatTemp <- surfWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% #Copies of water data
    left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1)) %>% mutate(predSet=makeTF(.,1-use))
    
  bWatTemp <- bottomWDat %>% unite(ID,c('YEID','doy'),sep='-') %>% 
    left_join(sDatTemp,by='ID') %>% filter(!is.na(PC1)) %>% mutate(predSet=makeTF(.,1-use))
  
  sWatPredSet <- sWatTemp %>% filter(predSet) #Take only samples from prediction set
  bWatPredSet <- bWatTemp %>% filter(predSet) 
  
  sWatTemp <- sWatTemp %>% filter(!predSet) #Remove prediction set
  bWatTemp <- bWatTemp %>% filter(!predSet)
  
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
  
  #Predict on withheld dataset, get difference
  sModPreds <- sWatPredSet %>% mutate(pred=predict(sMod,newdata=sWatPredSet),diff=pred-DO)
  bModPreds <- bWatPredSet %>% mutate(pred=predict(bMod,newdata=bWatPredSet),diff=pred-DO)
  
  return(list(surfaceDiff=sModPreds$diff,bottomDiff=bModPreds$diff))
}

#Function to run this in parallel, return as a dataframe
cvPredErrs <- function(i,inter=FALSE){
  source('helperFunctions.R')
  lapply(replicate(n=1000,expr=fitLagModsCV(i,interaction = inter),simplify = FALSE),function(x){
    ret <- data.frame(rmseSurface=rmse(x$surfaceDiff),rmseBottom=rmse(x$bottomDiff),
                      maeSurface=mae(x$surfaceDiff),maeBottom=mae(x$bottomDiff))
    return(ret)
  })
}

#Get prediction errors for surface/bottom, using RMSE and MAE

# library(parallel)
# cluster <- makeCluster(15)
# clusterExport(cluster,c('fitLagModsCV','sDat','surfWDat','bottomWDat'))
# cvPredList <- parLapply(cl=cluster,lags,cvPredErrs) #No interactions
# cvPredList2 <- parLapply(cl=cluster,lags,cvPredErrs,inter=TRUE) #Interactions
# stopCluster(cluster)
# cvPredList <- lapply(cvPredList,function(x) do.call('rbind',x))
# cvPredList2 <- lapply(cvPredList2,function(x) do.call('rbind',x))
# 
# cvPredList <- cvPredList %>% bind_rows(.id='lag') %>% mutate(lag=as.numeric(lag)) %>% pivot_longer(-lag) %>% 
#   mutate(depth=ifelse(grepl('Surface',name),'surface','bottom')) %>% 
#   mutate(errType=ifelse(grepl('rmse',name),'rmse','mae')) %>% select(-name) %>% mutate(modType='Simple')
# 
# cvPredList2 <- cvPredList2 %>% bind_rows(.id='lag') %>% mutate(lag=as.numeric(lag)) %>% pivot_longer(-lag) %>% 
#   mutate(depth=ifelse(grepl('Surface',name),'surface','bottom')) %>% 
#   mutate(errType=ifelse(grepl('rmse',name),'rmse','mae')) %>% select(-name) %>% mutate(modType='Interaction')
# 
# save(cvPredList,cvPredList2,file='./data/cvPredLists.Rdata')

load('./data/cvPredLists.Rdata')

bind_rows(cvPredList,cvPredList2) %>% filter(depth=='bottom') %>% 
  mutate(errType=factor(errType,labels=c('MAE','RMSE'))) %>% 
  ggplot(aes(x=lag,y=value,col=modType))+
  geom_point(alpha=0.1,position=position_dodge(width=0.5))+
  facet_wrap(~errType)+
  labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged linear Model',col='Model Type')+
  coord_cartesian(ylim=c(NA,5))

cvPredList %>% filter(depth=='bottom') %>% group_by(lag,errType) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
  mutate(errType=factor(errType,labels=c('MAE','RMSE'))) %>% 
  ggplot(aes(x=lag,y=med))+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.3)+
  geom_line()+facet_wrap(~errType)+
  labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged linear Model')

bind_rows(cvPredList,cvPredList2) %>% filter(depth=='bottom') %>% 
  mutate(errType=factor(errType,labels=c('MAE','RMSE'))) %>% 
  group_by(lag,errType,modType) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
  ggplot(aes(x=lag,y=med))+geom_ribbon(aes(ymax=max,ymin=min,fill=modType),alpha=0.3)+
  geom_line(aes(col=modType))+facet_wrap(~errType)+
  labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged linear Model',col='Model Type',fill='Model Type')+
  coord_cartesian(ylim=c(NA,10))
  
  

# Fit FR model of DO to gap-filled PCs ------------------------------------

NdayLag <- 30 #30 days in past
NdayForward <- 0 #0 days in future
dayLags <- -NdayForward:NdayLag

#Matrices to store PCA predictions for past 0:30 days
predMat <- matrix(NA,nrow=nrow(bottomWDat),ncol=length(dayLags),
                  dimnames=list(bottomWDat$YEID,gsub('-','m',paste0('lag',dayLags))))
pcaMatList <- list(PCA1=predMat,PCA2=predMat,PCA3=predMat,PCA4=predMat,PCA5=predMat)

for(p in 1:5){ #PCA dimensions
  for(i in 1:nrow(bottomWDat)){ #For each bottom water measurement
    getDays <- bottomWDat$doy[i]:bottomWDat$doy[i]-dayLags #Which days are 0-30 days behind the measurement?
    pcaMatList[[p]][i,] <- sDat %>% filter(sDat$YEID == bottomWDat$YEID[i] & sDat$doy %in% getDays) %>% pull(paste0('PC',p)) 
  }
}

#Data for functional regression
fdat <- list(DO_bottom=bottomWDat$DO,
             # DO_surf=surfWDat$DO,
             dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
             pcaMat1=pcaMatList$PCA1,pcaMat2=pcaMatList$PCA2,
             pcaMat3=pcaMatList$PCA3,pcaMat4=pcaMatList$PCA4,
             pcaMat5=pcaMatList$PCA5,doy=bottomWDat$doy,
             sE=bottomWDat$sE,sN=bottomWDat$sN,
             maxDepth=bottomWDat$maxDepth)

basisType <- 'cr' #Cubic regression splines
#Fit FDA models 
bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType)+s(dayMat,by=pcaMat5,bs=basisType), 
               data=fdat) #Bottom water

summary(bWatMod) #R-squared of about 0.50
par(mfrow=c(2,2)); gam.check(bWatMod); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
plot(bWatMod,scheme=1,pages=1)

#Use smoothPred to get FR plots from each smoother
p1 <- lapply(1:5,function(i){
  d <- expand.grid(dayMat=0:30,p=1) #Dataframe
  names(d)[2] <- paste0('pcaMat',i) #Change name of by variable
  smoothPred(m=bWatMod,dat=d,whichSmooth=i) 
}) %>% set_names(paste0('PCA',1:5)) %>% bind_rows(.id='PC') %>% 
  select(-contains('pcaMat')) %>% 
  ggplot(aes(x=dayMat))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred))+facet_wrap(~PC)+geom_hline(yintercept=0,col='red',linetype='dashed')+
  labs(x='Day (lag)',y='Effect')

p2 <- data.frame(pred=predict(bWatMod),actual=fdat$DO_bottom) %>% 
  ggplot()+geom_point(aes(x=pred,y=actual))+
  geom_abline(intercept = 0, slope = 1)+
  labs(x='Predicted Bottom DO',y='Actual Bottom DO')

(p <- ggarrange(p1,p2,ncol=2))
ggsave('./figures/frPCA_gapfill.png',p,width=10,height=5)
ggsave('./figures/frPCA_gapfill2.png',p1,width=10,height=5)

# Cross-validation (use 70%, predict on 30%)

#Function to fit models on random 70% of data, return prediction differences for remaining 30%
#uses fDat from above, and samples within it
fitFRmodCV <- function(i,dat=fdat,use=0.7){ #i doesn't do anything, just for parallel processing

  require(mgcv)
  source('helperFunctions.R')
  
  #Copy of spectral data
  predSet <- makeTF(data.frame(dat$DO_bottom),(1-use))
  
  getRows <- function(x,choose) if(any(class(x)=='matrix')) return(x[choose,]) else return(x[choose])
  
  predDat <- lapply(fdat,getRows,choose=predSet) #Data to predict on
  dat <- lapply(fdat,getRows,choose=!predSet) #Data to fit with
  
  basisType <- 'cr' #Cubic regression splines
  #Fit FDA models 
  bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                   s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType)+s(dayMat,by=pcaMat5,bs=basisType), 
                 data=dat) #Bottom water
  
  #Predict on withheld dataset, get difference
  bModPreds <- predict(bWatMod,newdata = predDat)-predDat$DO_bottom
  
  ret <- c(rmse(bModPreds),mae(bModPreds))
  names(ret) <- c('RMSE','MAE')
  return(ret)
}

library(parallel)
cluster <- makeCluster(15)
clusterExport(cluster,c('fitFRmodCV','fdat'))
cvPredList3 <- parLapply(cl=cluster,1:1000,fitFRmodCV)  #Takes only a few seconds to run
stopCluster(cluster)

pivot_longer(bind_rows(cvPredList3),RMSE:MAE) %>% group_by(name) %>% 
  # summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
  ggplot(aes(x=value))+geom_histogram()+facet_wrap(~name)+
  labs(x='Out-of-Sample Error',title='Bottom DO - Functional Regression Model')

#Compare models --------------------------

#RMSE
min(sapply(modList1,function(i) rmse(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) rmse(i$bottom))) #Lagged linear regression - PCA with interactions
rmse(bWatMod) #Functional regression - PCA

#Which days do these occur on?
lags[which.min(sapply(modList1,function(i) rmse(i$bottom)))] #7-day lag
lags[which.min(sapply(modList2,function(i) rmse(i$bottom)))] #0-day lag

#MAE
min(sapply(modList1,function(i) mae(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) mae(i$bottom))) #Lagged linear regression - PCA with interactions
mae(bWatMod) #Functional regression - PCA

#R2
max(sapply(modList1,function(i) getR2(i$bottom))) #Lagged linear regression - PCA
max(sapply(modList2,function(i) getR2(i$bottom)),na.rm=TRUE) #Lagged linear regression - PCA
summary(bWatMod)$r.sq #Functional regression - PCA

#df
sapply(modList1,function(i) getDF(i$bottom))[which.min(sapply(modList1,function(i) rmse(i$bottom)))]
sapply(modList2,function(i) getDF(i$bottom))[which.min(sapply(modList2,function(i) rmse(i$bottom)))]
bWatMod$df.residual


#Compare lagged linear to FR model
frStats <- pivot_longer(bind_rows(cvPredList3),RMSE:MAE) %>% group_by(name) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
  rename(errType=name) %>% mutate(errType=factor(errType,labels=c('MAE','RMSE')))

cvPredList %>% filter(depth=='bottom') %>% group_by(lag,errType) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
  mutate(errType=factor(errType,labels=c('MAE','RMSE'))) %>% 
  ggplot(aes(x=lag,y=med))+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.3)+
  geom_line()+
  geom_hline(data=frStats,aes(yintercept = med),col='red')+
  geom_hline(data=frStats,aes(yintercept = max),col='red',linetype='dashed')+
  geom_hline(data=frStats,aes(yintercept = min),col='red',linetype='dashed')+
  facet_wrap(~errType)+
  labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged Linear vs Functional Regression')
