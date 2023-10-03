## This script contains code to fit all the relative risk models. Note that fitting can take a long time, especially
## model L1. Also, the fitting function should be run a few times to improve model fit (see code below)
## The output is saved to separate R objects and used for the computation of absolute risk in a different script
## The example data was obtained by resampling with replacement and distorting ages and years at HL diagnosis for privacy reasons

# Loading the necessary scripts
source("RRfuns.R")


# Loading and processing data

load(file="20230904_exampledata.RData")



exampledata$agehlcat <- cut(exampledata$agehl, breaks=c(7,20,25,30,45), right=F)

exampledata$txyearcat <- cut(exampledata$txyear, breaks=c(1960, 1981,2005), right=F)

exampledata$ageyearstrat <- paste0(exampledata$agehlcat, exampledata$txyearcat)

# The analytic cohort
cohdat <- exampledata[exampledata$use_for_cohort_analysis==1,]

# NCC data
ccdat <- exampledata[exampledata$use_for_cc_analysis==1,]



# Computation of weights

ccdat$setsize <- sapply(ccdat$setnr, function(setn) sum(ccdat$setnr==setn))


ccdat$nriskstrat <- NA
ccdat$samwgtstrat<- NA
ccdat$setsizestrat  <- NA
ccdat$samwgt <- NA


for(mystrat in unique(ccdat$ageyearstrat)){
  
  ccdat[ccdat$ageyearstrat==mystrat,]$setsizestrat <- sapply(1:nrow(ccdat[ccdat$ageyearstrat==mystrat,]), function(i) sum(ccdat[ccdat$ageyearstrat==mystrat,]$setnr==ccdat[ccdat$ageyearstrat==mystrat,]$setnr[i] & ccdat[ccdat$ageyearstrat==mystrat,]$use_for_cohort_analysis==1))
  
  
  ccdat[ccdat$ageyearstrat==mystrat,]$nriskstrat <- sapply(1:nrow(ccdat[ccdat$ageyearstrat==mystrat,]), function(i) sum(cohdat[cohdat$ageyearstrat==mystrat,]$entry<=ccdat[ccdat$ageyearstrat==mystrat,]$exit[i] & cohdat[cohdat$ageyearstrat==mystrat,]$exit>=ccdat[ccdat$ageyearstrat==mystrat,]$exit[i]))
  
  
  ccdat[ccdat$ageyearstrat==mystrat,]$samwgtstrat <- sapply(1:nrow(ccdat[ccdat$ageyearstrat==mystrat,]), function(i) 1/(1-prod(1-(ccdat[ccdat$ageyearstrat==mystrat,][ccdat[ccdat$ageyearstrat==mystrat,]$exit > ccdat[ccdat$ageyearstrat==mystrat,]$entry[i]-3 & ccdat[ccdat$ageyearstrat==mystrat,]$exit <= ccdat[ccdat$ageyearstrat==mystrat,]$exit[i]+3 & ccdat[ccdat$ageyearstrat==mystrat,]$cc==1 & ccdat[ccdat$ageyearstrat==mystrat,]$use_for_cohort_analysis==1,]$setsizestrat-1)/(ccdat[ccdat$ageyearstrat==mystrat,][ccdat[ccdat$ageyearstrat==mystrat,]$exit > ccdat[ccdat$ageyearstrat==mystrat,]$entry[i]-3 & ccdat[ccdat$ageyearstrat==mystrat,]$exit <= ccdat[ccdat$ageyearstrat==mystrat,]$exit[i]+3 & ccdat[ccdat$ageyearstrat==mystrat,]$cc==1 & ccdat[ccdat$ageyearstrat==mystrat,]$use_for_cohort_analysis==1,]$nriskstrat-1))))
  
  ccdat[ccdat$ageyearstrat==mystrat & ccdat$use_for_cohort_analysis==0,]$samwgtstrat <- NA
  
  
  ccdat[ccdat$ageyearstrat==mystrat & ccdat$cc==1 & ccdat$use_for_cohort_analysis==1,]$samwgtstrat <- sum(cohdat$ageyearstrat==mystrat & cohdat$outcome=="BC")/sum(ccdat[ccdat$ageyearstrat==mystrat & ccdat$cc==1 & ccdat$use_for_cohort_analysis==1,]$cc)
  
  
  ccdat[ccdat$ageyearstrat==mystrat & ccdat$use_for_cohort_analysis==1,]$samwgt <- ccdat[ccdat$ageyearstrat==mystrat & ccdat$use_for_cohort_analysis==1,]$samwgtstrat
  
  
  ccdat[ccdat$ageyearstrat==mystrat & ccdat$cc==0 & ccdat$use_for_cohort_analysis==1,]$samwgt <- ccdat[ccdat$ageyearstrat==mystrat & ccdat$cc==0 & ccdat$use_for_cohort_analysis==1,]$samwgt*(sum(cohdat$ageyearstrat==mystrat & cohdat$outcome!="BC"))/sum(ccdat[ccdat$ageyearstrat==mystrat & ccdat$cc==0 & ccdat$use_for_cohort_analysis==1,]$samwgt)
  
  
}



fitdat <- ccdat


fitdat$agehlcat <- cut(fitdat$agehl, breaks=c(7,20,30,45), right=F)
fitdat$agehlcatnum <- as.numeric(fitdat$agehlcat)-1

currentdat <- cbind(fitdat[,c("setnr","cc","tumorQuadrant")], fitdat[,5:14])
currentdat <- cbind(currentdat,  model.matrix(~fam_bc_tot+Menopausal+Menopause_cat_num+AGEFLB+Parity_YN+agehlcatnum+txyearcat+use_for_cohort_analysis+samwgt, model.frame(~ ., fitdat, na.action=na.pass))[,-1])



#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Location-specific models




# L1

# Initial fit
fitL1 <- linearERRcomposite(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7, repar=TRUE,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="CCAL", weights=22)

# Euclidean norm of the score
sum(fitL1$MLE$score^2)


# Some more iterations to improve the fit
for(kk in 1:15){
  
  fitL1 <- linearERRcomposite(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7, repar=TRUE,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="CCAL", initpars=as.numeric(fitL1$MLE$fitobj$fit@fullcoef), weights=22)
  
  print(paste(sum(fitL1$MLE$score^2), fitL1$MLE$fitobj$fit@min, paste(as.numeric(round(fitL1$MLE$coef,3)), collapse=" ")))
}


save(fitL1, file="L1fit_example.RData")

#-------------

# L2

# Initial fit
fitL2 <- linearERRcompositeexp(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="CCAL", weights=22)

# Euclidean norm of the score
sum(fitL2$MLE$score^2)


# Some more iterations to improve the fit
for(kk in 1:15){
  
  fitL2 <- linearERRcompositeexp(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="CCAL", initpars=as.numeric(fitL2$MLE$fitobj$fit@fullcoef), weights=22)
  
  print(paste(sum(fitL2$MLE$score^2), fitL2$MLE$fitobj$fit@min, paste(as.numeric(round(fitL2$MLE$coef,3)), collapse=" ")))
}


save(fitL2, file="L2fit_example.RData")


#-------------

# L3

# Initial fit
fitL3 <- linearERRcompositecat(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="CCAL", weights=22)

# Euclidean norm of the score
sum(fitL3$MLE$score^2)


# Some more iterations to improve the fit
for(kk in 1:15){
  
  fitL3 <- linearERRcompositecat(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="CCAL", initpars=as.numeric(fitL3$MLE$fitobj$fit@fullcoef), weights=22)
  
  print(paste(sum(fitL3$MLE$score^2), fitL3$MLE$fitobj$fit@min, paste(as.numeric(round(fitL3$MLE$coef,3)), collapse=" ")))
}


save(fitL3, file="L3fit_example.RData")




#------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Mean dose models




# P1

currentdat[,4:13] <- ccdat$DMEAN

# Initial fit
fitP1 <- linearERRcomposite(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7, repar=TRUE,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="meandose", weights=22)

# Euclidean norm of the score
sum(fitP1$MLE$score^2)


# Some more iterations to improve the fit
for(kk in 1:15){
  
  fitP1 <- linearERRcomposite(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7, repar=TRUE,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="meandose", initpars=as.numeric(fitP1$MLE$fitobj$fit@fullcoef), weights=22)
  
  print(paste(sum(fitP1$MLE$score^2), fitP1$MLE$fitobj$fit@min, paste(as.numeric(round(fitP1$MLE$coef,3)), collapse=" ")))
}


save(fitP1, file="P1fit_example.RData")

#-------------

# P2
currentdat[,4:13] <- ccdat$DMEAN

# Initial fit
fitP2 <- linearERRcompositeexp(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="meandose", weights=22)

# Euclidean norm of the score
sum(fitP2$MLE$score^2)


# Some more iterations to improve the fit
for(kk in 1:15){
  
  fitP2 <- linearERRcompositeexp(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="meandose", initpars=as.numeric(fitP2$MLE$fitobj$fit@fullcoef), weights=22)
  
  print(paste(sum(fitP2$MLE$score^2), fitP2$MLE$fitobj$fit@min, paste(as.numeric(round(fitP2$MLE$coef,3)), collapse=" ")))
}


save(fitP2, file="P2fit_example.RData")


#-------------

# P3
currentdat[,4:13] <- ccdat$DMEAN

# Initial fit
fitP3 <- linearERRcompositecat(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="meandose", weights=22)

# Euclidean norm of the score
sum(fitP3$MLE$score^2)


# Some more iterations to improve the fit
for(kk in 1:15){
  
  fitP3 <- linearERRcompositecat(currentdat,4:13,set=1,status=2,loc=3,in_cohort=21,corrvars=14:20,corrvars_exclude_conditional=6:7,fitopt=list(maxit=100000, reltol=1e-12),ccmethod="meandose", initpars=as.numeric(fitP3$MLE$fitobj$fit@fullcoef), weights=22)
  
  print(paste(sum(fitP3$MLE$score^2), fitP3$MLE$fitobj$fit@min, paste(as.numeric(round(fitP3$MLE$coef,3)), collapse=" ")))
}


save(fitP3, file="P3fit_example.RData")




