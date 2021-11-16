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
load('./data/all2014_2.Rdata')

#GAM imputation ------------------------------------------

#Load smoothers
# load('./data/PCmods.RData') #Thin-plate splines
load('./data/PCmodsSoap.RData') #Soap film

#Get predictions of PCs at all locations through the entire season. If PC values missing, fill using PC model
sDat <- sDat %>% mutate(predPC1=predict(PCmod1,newdata=.),predPC2=predict(PCmod2,newdata=.),predPC3=predict(PCmod3,newdata=.)) %>% 
  mutate(predPC4=predict(PCmod4,newdata=.),predPC5=predict(PCmod5,newdata=.),predPC6=predict(PCmod6,newdata=.)) %>% 
  mutate(PC1=ifelse(gap,predPC1,PC1),PC2=ifelse(gap,predPC2,PC2),PC3=ifelse(gap,predPC3,PC3)) %>% 
  mutate(PC4=ifelse(gap,predPC4,PC4),PC5=ifelse(gap,predPC5,PC5),PC6=ifelse(gap,predPC6,PC6)) %>% 
  select(-predPC1:-predPC6) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j')) %>% #~10 seconds
  rename(imputed=gap)


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
    #Fit simple linear models use PC1:6
    sMod <- lm(DO~PC1+PC2+PC3+PC4+PC5+PC6,data=sWatTemp)
    bMod <- lm(DO~PC1+PC2+PC3+PC4+PC5+PC6,data=bWatTemp)
    
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

#Similar to models using raw data: no global minimum appears, but local min at ~12 days
#Simple model
(bestDay <- which.min(sapply(modList1,function(i) mae(i$bottom)))) #Minimum mae occurs on day 0 (no lag)
which.min(sapply(modList1,function(i) rmse(i$bottom))) #Minimum rmse occurs on day 0
#Interaction model
(bestDay <- which.min(sapply(modList2,function(i) mae(i$bottom)))) #Minimum mae occurs on day 3
which.min(sapply(modList2,function(i) rmse(i$bottom))) #Day 14

m1 <- modList1[[bestDay]]$bottom #Save model from that day
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
# cluster <- makeCluster(15)
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

bind_rows(cvPredList,cvPredList2) %>%  #Interaction model clearly has something wrong with it (probaby overfitting)
  ggplot(aes(x=lag,y=value,col=modType))+
  geom_point(alpha=0.1,position=position_dodge(width=0.5))+
  facet_wrap(~errType)+
  labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged linear Model',col='Model Type')+
  coord_cartesian(ylim=c(0,5))

cvPredList %>% filter(errType!='R2') %>% 
  group_by(lag,errType) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value),iqr=IQR(value)) %>% 
  ggplot(aes(x=lag,y=med))+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.3)+
  geom_line()+facet_wrap(~errType,ncol=1,scales='free_y')+
  labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged linear Model')

bind_rows(cvPredList,cvPredList2) %>% 
  group_by(lag,errType,modType) %>% 
  summarize(med=mean(value)) %>% ungroup() %>% 
  pivot_wider(names_from=errType,values_from = c(med)) %>% 
  group_by(modType) %>% 
  mutate(minMAE=MAE==min(MAE),minRMSE=RMSE==min(RMSE)) %>% 
  filter(minMAE|minRMSE) %>% select(-contains('min')) %>% 
  data.frame() %>% 
  relocate(lag,modType,RMSE,MAE,R2) 
  

# Fit FR model of DO to gap-filled PCs ------------------------------------

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
    pcaMatList[[p]][i,] <- sDat %>% filter(sDat$YEID == bottomWDat$YEID[i] & sDat$doy %in% getDays) %>% pull(paste0('PC',p)) 
  }
}

#Data for functional regression
fdat <- list(DO_bottom=bottomWDat$DO,
             # DO_surf=surfWDat$DO,
             dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
             pcaMat1=pcaMatList$PCA1,pcaMat2=pcaMatList$PCA2,
             pcaMat3=pcaMatList$PCA3,pcaMat4=pcaMatList$PCA4,
             pcaMat5=pcaMatList$PCA5,pcaMat6=pcaMatList$PCA6,
             doy=bottomWDat$doy,
             sE=bottomWDat$sE,sN=bottomWDat$sN,
             maxDepth=bottomWDat$maxDepth)

basisType <- 'cr' #Cubic regression splines
#Fit FDA models 
bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType)+s(dayMat,by=pcaMat5,bs=basisType)+
                 s(dayMat,by=pcaMat6,bs=basisType), 
               data=fdat) #Bottom water
save(bWatMod,file = './data/funRegMod.Rdata')

summary(bWatMod) #R-squared of about 0.61
par(mfrow=c(2,2)); gam.check(bWatMod); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
plot(bWatMod,scheme=1,pages=1)


#Use smoothPred to get FR plots from each smoother
pvals <- unname(round(summary(bWatMod)$s.table[,4],3))
pvals <- ifelse(pvals==0,'<0.001',paste0('=',as.character(pvals)))

p1 <- lapply(1:6,function(i){
  d <- expand.grid(dayMat=0:30,p=1) #Dataframe
  names(d)[2] <- paste0('pcaMat',i) #Change name of by variable
  smoothPred(m=bWatMod,dat=d,whichSmooth=i)}) %>% 
  set_names(paste0('PCA',1:6,' (p',pvals,')')) %>% 
  bind_rows(.id='PC') %>% 
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
  predDat$pred <- predict(bWatMod,newdata = predDat)
  predDat$diff <- predDat$pred-predDat$DO_bottom
  
  ret <- with(predDat,c(rmse(diff),mae(diff),summary(lm(DO_bottom~pred,data=predDat))$r.squared))
  names(ret) <- c('RMSE','MAE','R2')
  return(ret)
}

library(parallel)
cluster <- makeCluster(15)
clusterExport(cluster,c('fitFRmodCV','fdat'))
cvPredList3 <- parLapply(cl=cluster,1:1000,fitFRmodCV)  #Takes only a few seconds to run
stopCluster(cluster)

pivot_longer(bind_rows(cvPredList3),RMSE:R2) %>% group_by(name) %>% 
  summarize(mean=mean(value),med=median(value),max=max(value),min=min(value),iqr=IQR(value))

#What is the shortest lag time that we could use for FR and get similar results?

#Fit FDA models with different lags (10 days - 30 days)
lags <- 10:30

#Function to get error from model of lagged data. default = within-sample, 0<test<1 = out of sample
getLagErr <- function(l,f,test=0,N=30){
  f$dayMat <- f$dayMat[,1:l]
  f$pcaMat1 <- f$pcaMat1[,1:l]
  f$pcaMat2 <- f$pcaMat2[,1:l]
  f$pcaMat3 <- f$pcaMat3[,1:l]
  f$pcaMat4 <- f$pcaMat4[,1:l]
  f$pcaMat5 <- f$pcaMat5[,1:l]
  f$pcaMat6 <- f$pcaMat6[,1:l]
  
  basisType <- 'cr' #Cubic regression splines
  
  if(test==0){
    #Fit FDA models
    b <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
               s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType)+s(dayMat,by=pcaMat5,bs=basisType)+
               s(dayMat,by=pcaMat6,bs=basisType),
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
      b <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                 s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType)+s(dayMat,by=pcaMat5,bs=basisType)+
                 s(dayMat,by=pcaMat6,bs=basisType),
               data=traindat)
      
      #Get difference between predicted/actual from test data
      res <- testdat$DO_bottom-predict(b,newdata=testdat)
      ret <- c(rmse(res),mae(res),summary(b)$r.sq)
      names(ret) <- c('RMSE','MAE','R2')
      return(ret)
    }))
  }
}

bWatModList_train <- lapply(lags,getLagErr,f=fdat) %>% #Training data (within-sample) error
  bind_rows() %>% mutate(lag=lags) %>%
  pivot_longer(RMSE:R2) 

bWatModList_test <- lapply(lags,getLagErr,f=fdat,test=0.3,N=100) #Testing data (out of sample) error - takes a few mins

bWatModList_test <- bWatModList_test %>% lapply(.,data.frame) %>% set_names(nm=as.character(lags)) %>%
  bind_rows(.id = 'lag') %>%
  mutate(lag=as.numeric(lag)) %>%
  pivot_longer(RMSE:R2) %>%
  group_by(lag,name) %>% summarize(med=median(value),upr=quantile(value,0.9),lwr=quantile(value,0.1)) 

fdaErr <- left_join(bWatModList_test,bWatModList_train,by=c('lag','name')) %>% 
  rename('withinSampErr'='value')

fdaErr %>% ggplot(aes(x=lag,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+ #Testing
  geom_line(aes(y=withinSampErr),col='red')+
  facet_wrap(~name,ncol=1,scales='free_y')+
  labs(x='Maximum lag days',y='Model Accuracy')  

#Conclusion: Data from 30 days previous is better, but not a huge amount better than, data from a 10-day span

#Compare within-sample performance of models --------------------------

#RMSE
min(sapply(modList1,function(i) rmse(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) rmse(i$bottom))) #Lagged linear regression - PCA with interactions
rmse(bWatMod) #Functional regression - PCA

#Which days do these occur on?
lags[which.min(sapply(modList1,function(i) rmse(i$bottom)))] #12-day lag
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

# #Compare lagged linear to FR model
# frStats <- pivot_longer(bind_rows(cvPredList3),RMSE:MAE) %>% group_by(name) %>% 
#   summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
#   rename(errType=name) %>% mutate(errType=factor(errType,labels=c('MAE','RMSE')))
# 
# cvPredList %>% group_by(lag,errType) %>% 
#   summarize(mean=mean(value),med=median(value),max=max(value),min=min(value)) %>% 
#   ggplot(aes(x=lag,y=med))+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.3)+
#   geom_line()+
#   geom_hline(data=frStats,aes(yintercept = med),col='red')+
#   geom_hline(data=frStats,aes(yintercept = max),col='red',linetype='dashed')+
#   geom_hline(data=frStats,aes(yintercept = min),col='red',linetype='dashed')+
#   facet_wrap(~errType)+
#   labs(x='Time lag',y='Out-of-Sample Error',title='Bottom DO - Lagged Linear vs Functional Regression')

#Compare within-sample to out-of-sample error, using only simple lagged-linear model
llModWithin <- data.frame(lag=lags,MAE=sapply(modList1,function(i) mae(i$bottom)),
                         RMSE=sapply(modList1,function(i) rmse(i$bottom)),R2=sapply(modList1,function(i) getR2(i$bottom))) %>% 
  pivot_longer(cols=-lag,names_to='errType',values_to='withinSampErr')

llModErr <- cvPredList %>% mutate(lag=lag-1) %>% group_by(lag,errType) %>% 
  dplyr::summarize(med=median(value),max=max(value),min=min(value),iqr=IQR(value)) %>% 
  left_join(llModWithin,by=c('lag','errType'))

p1 <- ggplot(llModErr,aes(x=lag,y=med))+geom_ribbon(aes(ymax=max,ymin=min),alpha=0.3)+
  geom_line()+geom_line(aes(y=withinSampErr),col='red')+facet_wrap(~errType,ncol=1,scales='free_y')+
  labs(x='Time lag',y='Model Accuracy',title='Lagged linear Model')

# #Results from FDA model

# #Old method
# fdaWithin <- data.frame(MAE=mae(bWatMod),RMSE=rmse(mae(bWatMod)),R2=summary(bWatMod)$r.sq) %>% 
#   pivot_longer(cols=everything())
# fdaErr <- pivot_longer(bind_rows(cvPredList3),RMSE:R2) %>% group_by(name) %>% mutate(med=median(value)) 
# p2 <- fdaErr %>% 
#   ggplot(aes(x=value))+#geom_freqpoly()+
#   geom_histogram(fill='black',alpha=0.3,binwidth = 0.02)+
#   geom_vline(aes(xintercept=med),col='black')+
#   geom_vline(data=fdaWithin,aes(xintercept=value),col='red')+
#   facet_wrap(~name,ncol=1,scales='free_y')+
#   labs(x='Model Accuracy',title='Functional Data Analysis',y='Count')+
#   coord_flip()

p2 <- fdaErr %>% ggplot(aes(x=lag,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+ #Testing
  geom_line(aes(y=withinErr),col='red')+
  facet_wrap(~name,ncol=1,scales='free_y')+
  labs(x='Maximum lag days',y='Model Accuracy',title='Functional Data Analysis')  

(p <- ggarrange(p1,p2,nrow=1))
ggsave('./figures/outOfSamp_comparison.png',p,width=8,height=8)

NOTE <- 'Data are quite different from each other, so had to save as separate objects. llModErr = Out-of-sample and within-sample error for lagged linear models. fdaErr = Out-of-sample and within-sample error for 100 random draws of FDA model.'

save(llModErr,fdaErr,NOTE,file = './data/errDat.Rdata')
