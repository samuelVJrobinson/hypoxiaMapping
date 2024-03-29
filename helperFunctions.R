# Helper functions -----------

midcut <- function(dat,b,...){ #Function to cut data and assign label as midpoint (not boundaries)
  x <- cut(dat,breaks=b,...) # cut the data into bins...
  if(length(b)==1){ #If single break value
    s <- seq(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE),length.out=b+1) #Sequence of boundary points
  } else { #If multiple break values
    s <- b; b <- length(s)
  }
  levels(x) <- round(rowMeans(cbind(s,lag(s)))[2:c(b+1)],1) #Average with lag, then assign as new levels
  return(x)
}

dateCut <- function(d,b,include.lowest=TRUE){ #Cut date object d into b discrete chunks, and return as date object
  #include.lowest = should date on left boundary be included in category 1 (otherwise returns NA)
  if(class(b)=='character' & length(b)>1){ 
    b <- as.numeric(as.Date(b)) 
  } else if(!(length(b)==1 & class(b)=='numeric')) {
    stop('b must be single number, or character vector of date breaks in YYYY-MM-DD format')
  }
  cdates <- midcut(as.numeric(d),b,include.lowest=include.lowest) #Convert date to numeric, and chop into pieces
  cdates <- round(as.numeric(as.character(cdates))) #Turn into integers (using midpoint between breaks)
  cdates <- as.Date(cdates,origin="1970-01-01") #Turn into dates
  return(cdates)
}

#Turn sf geometry (x/y coordinates) into columns
#   removeGeom: should geometry be dropped?
#   epsgCode: should CRS be transformed (useful for lat/lon -> UTM)
geom2cols <- function(d,x=lon,y=lat,removeGeom=TRUE,epsg=NA){
  require(sf); require(dplyr)
  if(!('sf' %in% class(d))) stop('Data must be an sf object')
  if(is_grouped_df(d)){ #If data are grouped
    warning('Data are grouped. Will ungroup data before converting geometry column.')
    d <- ungroup(d)
  }
  if(!is.na(epsg)) d <- st_transform(d,epsg) #Transform to new CRS
  d <- d %>% mutate({{x}}:=st_coordinates(.)[,1],{{y}}:=st_coordinates(.)[,2]) #Make new columns from coordinates
  if(removeGeom) d <- st_drop_geometry(d) #Drop geometry
  return(d)
}


#Basic bivariate correlation plots (used in pairs() function)
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r),col=ifelse(r>0,'black','red'))
}

#Extract data for partial effects plots of smoothing terms
# dat = dataframe of predictor data (think "newdata" from predict.gam) + column of ones with title of "by" matrix (if doing functional regression)
#   * must have the same name as predictors
# m = GAM model
# whichSmooth = smoothing terms to use; all terms are aggregated (useful for interaction plots)
#   if numeric (vector/scalar), gets smooths in order stored in gam model. If character (vector/scalar), matches label exactly
# ci = multiplicitive factor for SE bounds (default = 1.96)
smoothPred <- function(dat,m,whichSmooth,ci=1.96){ 
  if(!(class(whichSmooth) %in% c('character','numeric','integer'))){
    stop('whichSmooth must be character, numeric, or integer')
  }
  
  if(is.character(whichSmooth)){ #Converts whichSmooth to numeric index, if in character form
    whichSmooth <- which(whichSmooth %in% sapply(m$smooth,function(x) x$label))
  }
  
  #Are "by" variables missing?
  byVars <- sapply(whichSmooth,function(x) m$smooth[[x]]$by)
  missByVar <- any(!byVars[byVars!='NA'] %in% names(dat))
  if(missByVar) stop(paste0('by variable(s) not specified. Possible names: ',byVars[byVars!='NA']))
  
  #Predictor matrices
  predMat <- lapply(whichSmooth,function(i) PredictMat(m$smooth[[i]],data=dat)) #Get predictor matrices from each smoother, using dat
  predMat <- do.call('cbind',predMat) #Amalgamate into single matrix
  #Coefficients
  coefRange <- do.call('c',lapply(m$smooth[whichSmooth],function(x) x$first.para:x$last.para)) #Get coefficients to use
  coefs <- coef(m)[coefRange] #Extract coefficient values
  #Predicted value - predictor matrix X coefficients
  dat$pred <- predMat %*% coefs
  #SE - swiped from plot.gam code
  dat$se <- sqrt(pmax(0,rowSums(predMat %*% m$Vp[coefRange,coefRange] * predMat)))
  #Confidence intervals
  dat$upr <- dat$pred + dat$se*ci; dat$lwr <- dat$pred - dat$se*ci 
  return(dat) #Return entire dataframe
}

#Changes negative values to scaleVal proportion of lowest positive value
fixNeg <- function(x,scaleVal=0.5) ifelse(x<0,scaleVal*min(x[x>0],na.rm=TRUE),x)
#Adds minimum value + minimum values * scaleVal 
fixNeg2 <- function(x,scaleVal=(1+1e-3)) x+abs(min(x,na.rm=TRUE))*scaleVal
#Linear rescaling between lwr and upr
rescale <- function(x,lwr=NA,upr=NA){
  if(is.na(lwr)&is.na(upr)) return(x)
  if(is.na(lwr)) lwr <- min(x,na.rm=TRUE)
  if(is.na(upr)) upr <- max(x,na.rm=TRUE)
  oldlwr <- min(x,na.rm=TRUE)
  oldupr <- max(x,na.rm=TRUE)
  oldcenter <- mean(c(oldlwr,oldupr))
  oldrng <- oldupr-oldlwr
  x <- (x-oldcenter)/oldrng
  center <- mean(c(lwr,upr)) #New midpoint
  rng <- upr-lwr #New range
  x <- (x*rng)+center #New data
  return(x)
}

#Mean absolute error
mae <- function(mod,backTrans=NULL){
  if(any(c('lm','gam') %in% class(mod))){
    res <- resid(mod)
  } else if(any(c('numeric','array') %in% class(mod))){
    res <- mod
  } else {
    return(NA)
  }
  if(!is.null(backTrans)){ #Transform before summation
    res <- eval(call(backTrans,res))
  }
  mean(abs(res))
}

#Root mean square error
rmse <- function(mod,backTrans=NULL){
  if(any(c('lm','gam') %in% class(mod))){
    res <- resid(mod)
  } else if(any(c('numeric','array') %in% class(mod))){
    res <- mod
  } else {
    return(NA)
  }
  if(!is.null(backTrans)){ #Transform before summation
    res <- eval(call(backTrans,res))
  }
  sqrt(mean(res^2))
} 

#R2 from lm objects
getR2 <- function(mod) if(class(mod)!='lm') return(NA) else summary(mod)$r.squared

#DF
getDF <- function(mod) if(class(mod)!='lm') return(NA) else mod$df.residual 
  
#Returns vector of random T/F of same length as dataframe d, where p = prob of T
makeTF <- function(d,p){ 
  N <- nrow(d)
  numT <- round(N*p)
  numF <- N-numT
  tf <- sample(c(rep(T,numT),rep(F,numF)))
  return(tf)
}  

#Get predictions from reduced model storage (list with coefs, vcv, and smooths)
predModList <- function(mList,newdat,returnSE=FALSE){
  #Predictor matrix
  predMat <- lapply(mList$smooths,function(x){
    PredictMat(x,data=newdat)}) %>% do.call('cbind',.)  
  predMat <- cbind(rep(1,nrow(predMat)),predMat) #Intercept column
  ret <- as.vector(predMat %*% mList$coefs) #Predictions
  if(returnSE){ #Standard error of prediction
    se <- sqrt(pmax(0,rowSums(predMat %*% mList$vcv * predMat)))  
    ret <- cbind(ret,se)
    colnames(ret) <- c('pred','se')
  }
  return(ret)
}
