#SECOND APPROACH TO ANALYSIS: MODEL BOTTOM DO BASED ON RAW SATELLITE CHANNELS
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

#Load data from saved file
load('./data/all2014.Rdata')

# Lagged linear model: PCA ----------------------------------------------------

lags <- 0:30 #0 to 30 day lags
fitLagMods <- function(i,interaction=FALSE){
  #Copy of spectral data
  sDatTemp <- st_drop_geometry(sDat) %>% 
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

# fitLagMods(7,TRUE)$surface
# debugonce(fitLagMods)

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

#Get plots of MSE and R-squared
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
p <- ggarrange(p1,p2,p3,ncol=1,common.legend=TRUE,legend='bottom') 
ggsave('./figures/lagPCAmod_raw.png',p,width=8,height=8)

# Functional regression using PCs -------------------------------------------------------------

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
    pcaMatList[[p]][i,] <- sDat %>% filter(sDat$YEID == bottomWDat$YEID[i] & sDat$doy %in% getDays) %>% pull(paste0('PC',p)) 
  }
}

#Data for functional regression
fdat <- list(DO_bottom=bottomWDat$DO,DO_surf=surfWDat$DO,
              dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
              pcaMat1=pcaMatList$PCA1,pcaMat2=pcaMatList$PCA2,
              pcaMat3=pcaMatList$PCA3,pcaMat4=pcaMatList$PCA4,
              doy=bottomWDat$doy,sE=bottomWDat$sE,sN=bottomWDat$sN,
              maxDepth=bottomWDat$maxDepth)

#Problem: lots of NAs. Some kind of imputation is necessary
plot(0:30,c(sum(!is.na(fdat$pcaMat1[,1])),
  sapply(2:ncol(fdat$pcaMat1),function(x){
  sum(!is.na(apply(fdat$pcaMat1[,c(1:x)],1,sum)))
})),ylab='Number of non-NA measurements',xlab='Day Lag',pch=19)



basisType <- 'cr' #Thin-plate regression splines with extra shrinkage. Cubic splines have higher R2 but have very strange shapes
#Fit FDA models 
bWatMod <- gam(DO_bottom ~ s(dayMat,by=pcaMat1,bs=basisType)+s(dayMat,by=pcaMat2,bs=basisType)+
                  s(dayMat,by=pcaMat3,bs=basisType)+s(dayMat,by=pcaMat4,bs=basisType), 
                data=fdat) #Bottom water

#Summary: FR isn't possible without some kind of missing data imputation, as there is too much missing data

# Compare models ---------------------------

#RMSE
min(sapply(modList1,function(i) rmse(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) rmse(i$bottom)),na.rm=TRUE) #Lagged linear regression - PCA with interactions

#Which days do these occur on?
lags[which.min(sapply(modList1,function(i) rmse(i$bottom)))] #11-day gap
lags[which.min(sapply(modList2,function(i) rmse(i$bottom)))] #13-day gap

#MAE
min(sapply(modList1,function(i) mae(i$bottom))) #Lagged linear regression - PCA
min(sapply(modList2,function(i) mae(i$bottom)),na.rm=TRUE) #Lagged linear regression - PCA with interactions

#R2
max(sapply(modList1,function(i) getR2(i$bottom))) #Lagged linear regression - PCA
max(sapply(modList2,function(i) getR2(i$bottom)),na.rm=TRUE) #Lagged linear regression - PCA

#df
sapply(modList1,function(i) getDF(i$bottom))[which.min(sapply(modList1,function(i) rmse(i$bottom)))]
sapply(modList2,function(i) getDF(i$bottom))[which.min(sapply(modList2,function(i) rmse(i$bottom)))]

