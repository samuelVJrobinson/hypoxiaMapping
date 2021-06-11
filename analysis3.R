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

pcDat <- sDat %>% select(YEID,date_img,gap,E:geometry) %>% filter(!gap) %>% #Strips out missing data
  mutate(sE=sE+rnorm(n(),0,0.1),sN=sN+rnorm(n(),0,0.1)) #Add a bit of random noise to make sure that distances between points aren't 0

PCmod1 <-gam(PC1~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat)
PCmod2 <-gam(PC2~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=pcDat)
PCmod3 <-gam(PC3~te(sN,sE,doy,bs=c('tp','tp'),k=c(40,10),d=c(2,1)),data=pcDat) 
PCmod4 <-gam(PC4~te(sN,sE,doy,bs=c('tp','tp'),k=c(40,10),d=c(2,1)),data=pcDat) #Not that great
PCmod5 <-gam(PC5~te(sN,sE,doy,bs=c('tp','tp'),k=c(40,10),d=c(2,1)),data=pcDat) #Not that great

# par(mfrow=c(2,2)); 
# gam.check(PCmod1); abline(0,1,col='red'); 
# gam.check(PCmod2); abline(0,1,col='red'); 
# gam.check(PCmod3); abline(0,1,col='red'); 
# gam.check(PCmod4); abline(0,1,col='red'); #Not as good, but minimal variance, so this probably doesn't matter as much
# gam.check(PCmod5); abline(0,1,col='red'); 
# par(mfrow=c(1,1))

#See how predictions look on spatial grid - takes a few minutes
predPC <- with(pcDat,expand.grid(doy=seq(min(doy),max(doy),length.out=9),loc=unique(paste(sE,sN,sep='_')))) %>%
  separate(loc,c('sE','sN'),sep='_',convert=TRUE) %>%
  mutate(predPC1=predict(PCmod1,newdata=.),sePC1=predict(PCmod1,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(predPC2=predict(PCmod2,newdata=.),sePC2=predict(PCmod2,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(predPC3=predict(PCmod3,newdata=.),sePC3=predict(PCmod3,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(predPC4=predict(PCmod4,newdata=.),sePC4=predict(PCmod4,newdata=.,se.fit=TRUE)$se.fit) %>%
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

#Predicted values
p1 <- predPC %>% select(date,sE,sN,contains('pred')) %>%
  pivot_longer(contains('pred')) %>% mutate(name=gsub('pred','',name)) %>%
  ggplot()+geom_point(aes(x=sE,y=sN,col=value))+
  facet_grid(name~date)+
  scale_colour_distiller(type='div',palette = "Spectral",direction=-1)+
  labs(x='E',y='N',col='PC Value',title='Predicted value')+
  theme(legend.position='bottom')
# 
# #SE of prediction
# p2 <- predPC %>% select(date,sE,sN,contains('sePC')) %>% 
#   pivot_longer(contains('sePC')) %>% mutate(name=gsub('sePC','PC',name)) %>% 
#   ggplot()+geom_point(aes(x=sE,y=sN,col=value))+
#   facet_grid(name~date)+
#   scale_colour_distiller(type='div',palette = "Reds",direction=1)+
#   labs(x='E',y='N',col='SE of PC Value',title='Standard error of Prediction')+
#   theme(legend.position='bottom')
# 
# #Residual plots
# p3 <- pcDat %>% mutate(residPC1=resid(PCmod1),residPC2=resid(PCmod2),residPC3=resid(PCmod3),residPC4=resid(PCmod4)) %>% 
#   pivot_longer(contains('residPC')) %>% mutate(name=gsub('residPC','PC',name)) %>% 
#   mutate(date=dateCut(date_img,9)) %>% 
#   ggplot()+geom_point(aes(x=sE,y=sN,col=value,alpha=abs(value),size=abs(value)))+
#   facet_grid(name~date)+
#   scale_colour_distiller(type='div',palette = "Spectral",direction=1)+
#   scale_alpha_continuous(range=c(0.01,0.75))+
#   labs(x='E',y='N',col='Residual',title='Residual plot')+
#   theme(legend.position='bottom')

#Get predictions of PCs at all locations through the entire season
sDat <- sDat %>% mutate(predPC1=predict(PCmod1,newdata=.),predPC2=predict(PCmod2,newdata=.)) %>% 
  mutate(predPC3=predict(PCmod3,newdata=.),predPC4=predict(PCmod4,newdata=.),predPC5=predict(PCmod5,newdata=.)) %>% 
  mutate(PC1=ifelse(gap,predPC1,PC1),PC2=ifelse(gap,predPC2,PC2),PC3=ifelse(gap,predPC3,PC3),PC4=ifelse(gap,predPC4,PC4),PC5=ifelse(gap,predPC4,PC5)) %>% 
  select(-predPC1:-predPC5) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j')) 

#Fit model of DO to gap-filled PCs ----------------------------------

lags <- 0:30 #Try 0 to 30 day lags
modList <- lapply(lags,function(i){ #Lagged linear models of DO using PC1:4
  #Copy of spectral data
  sDatTemp <- st_drop_geometry(sDat) %>% 
    select(YEID,doy,contains('PC')) %>% 
    mutate(doy=doy+i) #Add to go back, subtract to go forward
  
  sWatTemp <- surfWDat %>% #Copies of water data
    left_join(sDatTemp,by=c('YEID','doy'))
  bWatTemp <- bottomWDat %>% 
    left_join(sDatTemp,by=c('YEID','doy'))
  
  #Fit simple linear models use PC1:4
  sMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=sWatTemp)
  bMod <- lm(DO~PC1+PC2+PC3+PC4+PC5,data=bWatTemp)
  return(list(surface=sMod,bottom=bMod))
})

modList2 <- lapply(lags,function(i){ #Lagged linear models of DO using PC1:4 (interactions)
  #Copy of spectral data
  sDatTemp <- st_drop_geometry(sDat) %>% 
    select(YEID,doy,contains('PC')) %>% 
    mutate(doy=doy+i) #Add to go back, subtract to go forward
  
  sWatTemp <- surfWDat %>% #Copies of water data
    left_join(sDatTemp,by=c('YEID','doy'))
  bWatTemp <- bottomWDat %>% 
    left_join(sDatTemp,by=c('YEID','doy'))
  
  #Interactions of PCs
  sMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=sWatTemp)
  bMod <- lm(DO~PC1*PC2*PC3*PC4*PC5,data=bWatTemp)
  return(list(surface=sMod,bottom=bMod))
})

#Get plots of MSE and R-squared
p1 <- data.frame(lag=lags,
                 surface1=sapply(modList,function(i) mean(abs(resid(i$surface)))),
                 bottom1=sapply(modList,function(i) mean(abs(resid(i$bottom)))),
                 surface2=sapply(modList2,function(i) mean(abs(resid(i$surface)))),
                 bottom2=sapply(modList2,function(i) mean(abs(resid(i$bottom))))) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='Mean Absolute Error',col='Model\nType')

p2 <- data.frame(lag=lags,
                 surface1=sapply(modList,function(i) summary(i$surface)$r.squared),
                 bottom1=sapply(modList,function(i) summary(i$bottom)$r.squared),
                 surface2=sapply(modList2,function(i) summary(i$surface)$r.squared),
                 bottom2=sapply(modList2,function(i) summary(i$bottom)$r.squared)) %>% 
  pivot_longer(surface1:bottom2) %>% 
  mutate(modType=ifelse(grepl('1',name),'Simple','Interaction'),name=gsub('(1|2)','',name)) %>% 
  ggplot()+geom_line(aes(x=lag,y=value,col=modType))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='R-squared',col='Model\nType')
ggarrange(p1,p2,ncol=1,common.legend=TRUE,legend='right') 

#Similar to models using raw data: no local minimum appears

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
fdat <- list(DO_bottom=bottomWDat$DO,DO_surf=surfWDat$DO,
             dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
             pcaMat1=pcaMatList$PCA1,pcaMat2=pcaMatList$PCA2,
             pcaMat3=pcaMatList$PCA3,pcaMat4=pcaMatList$PCA4,
             pcaMat5=pcaMatList$PCA5,doy=bottomWDat$doy,
             sE=bottomWDat$sE,sN=bottomWDat$sN,
             maxDepth=bottomWDat$maxDepth)

basisType <- 'ts' #Thin-plate regression splines with extra shrinkage. Cubic splines have higher R2 but have very strange shapes
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

ggarrange(p1,p2,ncol=2)


