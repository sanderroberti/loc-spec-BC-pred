library(survival)
library(splines)


## This script illustrates how absolute risks are computed for patients in the case-control study
## Also, this script includes the fitting of the competing risk model
## The example data was obtained by resampling with replacement and distorting ages and years at HL diagnosis for privacy reasons


# Load one or multiple of the relative risk models
load(file="L1fit_example.RData")
load(file="L2fit_example.RData")
load(file="L3fit_example.RData")
load(file="P1fit_example.RData")
load(file="P2fit_example.RData")
load(file="P3fit_example.RData")

load(file="genpopinc.RData")
load(file="survinc.RData")

mymod <- "L1" # Model to use to predict risks. Options are L1, L2, L3, P1, P2, and P3


# Loading and preparing data

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




subdata <- ccdat
subdata$agehlcat <- cut(subdata$agehl, breaks=c(7,20,30,45), right=F)
subdata$agehlcatnum <- as.numeric(subdata$agehlcat)-1

currentdat <- cbind(subdata[,c("setnr","cc","tumorQuadrant")], subdata[,5:14])
currentdat <- cbind(currentdat,  model.matrix(~fam_bc_tot+Menopausal+Menopause_cat_num+AGEFLB+Parity_YN+agehlcatnum+txyearcat+use_for_cohort_analysis+samwgt, model.frame(~ ., subdata, na.action=na.pass))[,-1])

# Computing relative risk at NCC study index date for the purpose of estimating phi
if(mymod == "L1"){
  # L1
  RRmod <- fitL1
  rr_AR <- c(exp(fitL1$MLE$coef[7:13]%*%t(as.matrix(currentdat[,14:20]))))*
    rowSums((1+exp(fitL1$MLE$coef[2])*currentdat[,4:13])*
              matrix(rep(rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)), nrow(currentdat)), nrow=nrow(currentdat), byrow=TRUE))#/
  
} else if (mymod=="L2"){
  # L2
  RRmod <- fitL2
  rr_AR <- c(exp(fitL2$MLE$coef[7:13]%*%t(as.matrix(currentdat[,14:20]))))*
    rowSums(exp(fitL2$MLE$coef[2]*currentdat[,4:13])*
              matrix(rep(rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)), nrow(currentdat)), nrow=nrow(currentdat), byrow=TRUE))
  
} else if (mymod=="L3"){
  # L3
  RRmod <- fitL3
  rr_AR <- c(exp(fitL3$MLE$coef[8:14]%*%t(as.matrix(currentdat[,14:20]))))*
    rowSums(sapply(currentdat[,4:13], function(x) exp(model.matrix(~cut(x, breaks=c(0,10,25,500),right=FALSE)-1)[,-1]%*%fitL3$MLE$coef[2:3]))*
              matrix(rep(rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)/sum(rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)), nrow(currentdat)), nrow=nrow(currentdat), byrow=TRUE))#/
  
  
} else if (mymod=="P1"){
  # P1
  RRmod <- fitP1
  currentdat[,4:13] <- ccdat$DMEAN
  rr_AR <- c(exp(fitP1$MLE$coef[3:9]%*%t(as.matrix(currentdat[,14:20]))))*
    (1+exp(fitP1$MLE$coef[2])*rowMeans(currentdat[,4:13]))
  
} else if (mymod=="P2"){
  # P2
  RRmod <- fitP2
  currentdat[,4:13] <- ccdat$DMEAN
  rr_AR <- c(exp(fitP2$MLE$coef[3:9]%*%t(as.matrix(currentdat[,14:20]))))*
    exp(fitP2$MLE$coef[2]*rowMeans(currentdat[,4:13]))
  
} else if (mymod=="P3"){
  # P3
  RRmod <- fitP3
  currentdat[,4:13] <- ccdat$DMEAN
  rr_AR <- c(exp(fitP3$MLE$coef[4:10]%*%t(as.matrix(currentdat[,14:20]))))*
    exp(model.matrix(~cut(ccdat$DMEAN, breaks=c(0,10,25,500),right=FALSE))[,-1]%*%fitP3$MLE$coef[2:3])
  
}



phi_under55_gb <- 1-1/(sum(rr_AR[ccdat$exit<55 & ccdat$use_for_cohort_analysis==1]*subdata$samwgt[ccdat$exit<55 & ccdat$use_for_cohort_analysis==1])/sum(subdata$samwgt[ccdat$exit<55 & ccdat$use_for_cohort_analysis==1]))
phi_over55_gb <- 1-1/(sum(rr_AR[ccdat$exit>=55 & ccdat$use_for_cohort_analysis==1]*subdata$samwgt[ccdat$exit>=55 & ccdat$use_for_cohort_analysis==1])/sum(subdata$samwgt[ccdat$exit>=55 & ccdat$use_for_cohort_analysis==1]))

phi <- function(age){
  ifelse(age<55, phi_under55_gb,phi_over55_gb)
}

xs2 <- sort(unique(c(genpopinc$age,99, unique(ccdat[ccdat$Lft_eerste_kind_combi_cleaned>=min(genpopinc$age),]$Lft_eerste_kind_combi_cleaned), unique(ccdat[ccdat$Menop_age_combi>=min(genpopinc$age),]$Menop_age_combi) )))

genpopinc2 <- genpopinc$incidence[floor(xs2)-15]




if(mymod == "L1"){
  # L1
  rrmat <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitL1$MLE$coef[9] + 1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitL1$MLE$coef[8]+
          parityvar*fitL1$MLE$coef[10] + 1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitL1$MLE$coef[11])*
      c(exp(fitL1$MLE$coef[c(7,12,13)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      sum((1+exp(fitL1$MLE$coef[2])*currentdat[myrrrow,4:13])*
            rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)))#/
    
  }))
} else if (mymod=="L2"){
  # L2
  rrmat <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitL2$MLE$coef[9] + 1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitL2$MLE$coef[8]+
          parityvar*fitL2$MLE$coef[10] + 1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitL2$MLE$coef[11])*
      c(exp(fitL2$MLE$coef[c(7,12,13)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      sum((exp(fitL2$MLE$coef[2]*currentdat[myrrrow,4:13]))*
            rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)))#/
    
  }))
} else if (mymod=="L3"){
  # L3
  rrmat <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    exp(menopvar*fitL3$MLE$coef[10] + 1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitL3$MLE$coef[9]+
          parityvar*fitL3$MLE$coef[11] + 1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitL3$MLE$coef[12])*
      c(exp(fitL3$MLE$coef[c(8,13,14)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      sum((sapply(currentdat[myrrrow,4:13], function(x) exp(model.matrix(~cut(x, breaks=c(0,10,25,500),right=FALSE)-1)[,-1]%*%fitL3$MLE$coef[2:3])))*
            rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)/sum(rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)))#/
    
  }))
  
} else if (mymod=="P1"){
  # P1
  rrmat <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    exp(menopvar*fitP1$MLE$coef[5] + 1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitP1$MLE$coef[4]+
          parityvar*fitP1$MLE$coef[6] + 1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitP1$MLE$coef[7])*
      c(exp(fitP1$MLE$coef[c(3,8,9)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      (1+exp(fitP1$MLE$coef[2])*mean(as.numeric(currentdat[myrrrow,4:13])))
    
    
  }))
  
} else if (mymod=="P2"){
  # P2
  rrmat <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitP2$MLE$coef[5] + 1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitP2$MLE$coef[4]+
          parityvar*fitP2$MLE$coef[6] + 1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitP2$MLE$coef[7])*
      c(exp(fitP2$MLE$coef[c(3,8,9)]%*%t(as.matrix(currentdat[myrrrow,c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)")]))))*
      exp(fitP2$MLE$coef[2]*currentdat[myrrrow,4])
    
  }))
  
  
  
} else if (mymod=="P3"){
  # P3
  rrmat <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitP3$MLE$coef[6] + 1*(xs2>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitP3$MLE$coef[5]+
          parityvar*fitP3$MLE$coef[7] + 1*(xs2>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitP3$MLE$coef[8])*
      c(exp(fitP3$MLE$coef[c(4,9,10)]%*%t(as.matrix(currentdat[myrrrow,c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)")]))))*
      c(exp(model.matrix(~cut(ccdat[myrrrow,]$DMEAN, breaks=c(0,10,25,500),right=FALSE))[,-1]%*%fitP3$MLE$coef[2:3]))
    
  }))
}




rrmat <- rrmat[,-length(xs2)]


Tmat0 <- sapply(1:(length(xs2)-1), function(age){
  pmin(pmin(pmax(pmin(ccdat$exit,45)-xs2[age],0),xs2[age+1]-xs2[age]), pmin(pmax(xs2[age+1]-ccdat$entry,0),xs2[age+1]-xs2[age]))
})
Tmat1 <- sapply(1:(length(xs2)-1), function(age){
  pmin(pmin(pmax(ccdat$exit-xs2[age],0),xs2[age+1]-xs2[age]), pmin(pmax(xs2[age]-pmax(ccdat$entry,45)+1,0),xs2[age+1]-xs2[age]))
})


# Compute the integrated `baseline hazard' (without rho) for the computation of rho
ccdat$popinc_under45 <- c((rrmat*Tmat0)%*%(genpopinc2[-length(xs2)]*(1-phi(xs2[-length(xs2)]))))
ccdat$popinc_over45 <- c((rrmat*Tmat1)%*%(genpopinc2[-length(xs2)]*(1-phi(xs2[-length(xs2)]))))




survinc_nonzero <- survinc[survinc$hazard0>0,]
xs2surv <- sort(survinc_nonzero$age)




if(mymod == "L1"){
  # L1
  rrmatsurv <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2surv, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2surv, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitL1$MLE$coef[9] + 1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitL1$MLE$coef[8]+
          parityvar*fitL1$MLE$coef[10] + 1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitL1$MLE$coef[11])*
      c(exp(fitL1$MLE$coef[c(7,12,13)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      sum((1+exp(fitL1$MLE$coef[2])*currentdat[myrrrow,4:13])*
            rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)))#/
    
  }))
} else if (mymod=="L2"){
  # L2
  rrmatsurv <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2surv, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2surv, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitL2$MLE$coef[9] + 1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitL2$MLE$coef[8]+
          parityvar*fitL2$MLE$coef[10] + 1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitL2$MLE$coef[11])*
      c(exp(fitL2$MLE$coef[c(7,12,13)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      sum((exp(fitL2$MLE$coef[2]*currentdat[myrrrow,4:13]))*
            rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)))#/
    
  }))
} else if (mymod=="L3"){
  # L3
  rrmatsurv <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2surv, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2surv, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    exp(menopvar*fitL3$MLE$coef[10] + 1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitL3$MLE$coef[9]+
          parityvar*fitL3$MLE$coef[11] + 1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitL3$MLE$coef[12])*
      c(exp(fitL3$MLE$coef[c(8,13,14)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      sum((sapply(currentdat[myrrrow,4:13], function(x) exp(model.matrix(~cut(x, breaks=c(0,10,25,500),right=FALSE)-1)[,-1]%*%fitL3$MLE$coef[2:3])))*
            rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)/sum(rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)))#/
    
  }))
  
} else if (mymod=="P1"){
  # P1
  rrmatsurv <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2surv, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2surv, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    exp(menopvar*fitP1$MLE$coef[5] + 1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitP1$MLE$coef[4]+
          parityvar*fitP1$MLE$coef[6] + 1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitP1$MLE$coef[7])*
      c(exp(fitP1$MLE$coef[c(3,8,9)]%*%t(as.matrix(currentdat[myrrrow,c(c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)"))]))))*
      (1+exp(fitP1$MLE$coef[2])*mean(as.numeric(currentdat[myrrrow,4:13])))
    
    
  }))
  
} else if (mymod=="P2"){
  # P2
  rrmatsurv <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2surv, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2surv, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitP2$MLE$coef[5] + 1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitP2$MLE$coef[4]+
          parityvar*fitP2$MLE$coef[6] + 1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitP2$MLE$coef[7])*
      c(exp(fitP2$MLE$coef[c(3,8,9)]%*%t(as.matrix(currentdat[myrrrow,c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)")]))))*
      exp(fitP2$MLE$coef[2]*currentdat[myrrrow,4])
    
  }))
  
  
  
} else if (mymod=="P3"){
  # P3
  rrmatsurv <- do.call("rbind",lapply(1:nrow(ccdat), function(myrrrow){
    menopvar <- as.numeric(factor(cut(xs2surv, breaks=c(0, 29.99,39.99,49.99,100))))*1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)
    
    if(any(menopvar>0)) menopvar[menopvar>0] <- min(menopvar[menopvar>0])
    
    
    parityvar <- (as.numeric(factor(cut(xs2surv, breaks=c(0, 25,30,100), right=FALSE)))-1)*1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")
    
    if(any(parityvar>0)) parityvar[parityvar>0] <- min(parityvar[parityvar>0])
    
    
    
    exp(menopvar*fitP3$MLE$coef[6] + 1*(xs2surv>=ccdat[myrrrow,]$Menop_age_combi)*(ccdat[myrrrow,]$Menopausal==1)*fitP3$MLE$coef[5]+
          parityvar*fitP3$MLE$coef[7] + 1*(xs2surv>=ccdat[myrrrow,]$Lft_eerste_kind_combi_cleaned)*(ccdat[myrrrow,]$Parity_YN=="Ja")*fitP3$MLE$coef[8])*
      c(exp(fitP3$MLE$coef[c(4,9,10)]%*%t(as.matrix(currentdat[myrrrow,c("fam_bc_tot","agehlcatnum","txyearcat[1.98e+03,2e+03)")]))))*
      c(exp(model.matrix(~cut(ccdat[myrrrow,]$DMEAN, breaks=c(0,10,25,500),right=FALSE))[,-1]%*%fitP3$MLE$coef[2:3]))
    
  }))
}


ccdat$popinc_under45surv <-  sapply(1:nrow(rrmat), function(myrow){
  sum(rrmatsurv[myrow,xs2surv<45 & xs2surv<=ccdat[myrow,]$exit & xs2surv>ccdat[myrow,]$entry]*survinc_nonzero[survinc_nonzero$age<45 & survinc_nonzero$age<=ccdat[myrow,]$exit & survinc_nonzero$age>ccdat[myrow,]$entry,]$hazard0*(1-phi(xs2surv[xs2surv<45 & xs2surv<=ccdat[myrow,]$exit & xs2surv>ccdat[myrow,]$entry])))
})

ccdat$popinc_over45surv <-  sapply(1:nrow(rrmat), function(myrow){
  sum(rrmatsurv[myrow,xs2surv>=45 & xs2surv<=ccdat[myrow,]$exit & xs2surv>ccdat[myrow,]$entry]*survinc_nonzero[survinc_nonzero$age>=45 & survinc_nonzero$age<=ccdat[myrow,]$exit & survinc_nonzero$age>ccdat[myrow,]$entry,]$hazard0*(1-phi(xs2surv[xs2surv>=45 & xs2surv<=ccdat[myrow,]$exit & xs2surv>ccdat[myrow,]$entry])))
})




ccdat$agehl <- floor(ccdat$agehl)
ccdat$exit <- floor(ccdat$exit)
ccdat$entry <- floor(ccdat$entry)

cohdat$agehl <- floor(cohdat$agehl)
cohdat$exit <- floor(cohdat$exit)
cohdat$entry <- floor(cohdat$entry)



res2 <- ccdat[,c("entry","exit","setnr","cc","samwgt", "txyearcat","popinc_under45","popinc_over45", "popinc_under45surv","popinc_over45surv", "use_for_cohort_analysis") ]
res2 <- cbind(res2, RR=rr_AR)
res2 <- res2[order(res2$exit),]
res2$setsize <- sapply(res2$setnr, function(xset) sum(ccdat$setnr==xset))



rho_under45 <- sum((res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$exit < 45)*res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$samwgt)/sum(res2[ res2$use_for_cohort_analysis==1,]$popinc_under45*res2[ res2$use_for_cohort_analysis==1,]$samwgt)
rho_over45 <- sum((res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$exit >= 45)*res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$samwgt)/sum(res2[ res2$use_for_cohort_analysis==1,]$popinc_over45*res2[ res2$use_for_cohort_analysis==1,]$samwgt)


rhofun <- function(age){
  ifelse(age <45, rho_under45, rho_over45)
}


rho_under45surv <- sum((res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$exit < 45)*res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$samwgt)/sum(res2[ res2$use_for_cohort_analysis==1,]$popinc_under45surv*res2[ res2$use_for_cohort_analysis==1,]$samwgt)
rho_over45surv <- sum((res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$exit >= 45)*res2[res2$cc==1 & res2$use_for_cohort_analysis==1,]$samwgt)/sum(res2[ res2$use_for_cohort_analysis==1,]$popinc_over45surv*res2[ res2$use_for_cohort_analysis==1,]$samwgt)
rhofunsurv <- function(age){
  ifelse(age <45, rho_under45surv, rho_over45surv)
}




res <- ccdat[ccdat$cc==1,c("entry","exit","setnr","use_for_cohort_analysis")]
res <- cbind(res, RR=rr_AR[ccdat$cc==1])
res <- res[order(res$exit),]
res$natrisk <- sapply(res$exit, function(x) sum(cohdat[!is.na(cohdat$exit),]$entry <= x & cohdat[!is.na(cohdat$exit),]$exit >= x))



res2 <- res[res$use_for_cohort_analysis==1,]

res2$agecat <- cut(res2$exit, breaks=c(0,30,35,40,43,45,47,50,53,55,57,100), right=FALSE)
cohdat$agecat <- cut(cohdat$exit, breaks=c(0,30,35,40,43,45,47,50,53,55,57,100), right=FALSE)

res2$wgt <- sapply(res2$agecat, function(x) sum(cohdat[cohdat$agecat==x,]$outcome=="BC")/sum(res2$agecat==x))


myxs2 <- unique(res2$exit)
myys2 <- sapply(myxs2, function(x){
  subres <- res2[res2$exit==x,]
  sum(sapply(1:nrow(subres), function(xx){
    subres[xx,]$wgt/(mean(rr_AR[ccdat$setnr==subres[xx,]$setnr & ccdat$use_for_cohort_analysis==1])*subres[xx,]$natrisk)
  }))
})

lbbase <- data.frame(age=myxs2, h1_base=myys2) # Cohort-based baseline (Langholz-Borgan)



liubase <- data.frame(age=genpopinc$age, baseline=genpopinc$incidence*(1-phi(genpopinc$age))*rhofun(genpopinc$age)) # Cohort + general population baseline


survbase <- data.frame(age=survinc$age, baseline=survinc$hazard0*(1-phi(survinc$age))*rhofunsurv(survinc$age)) # Cohort + survivors in general population baseline


#-------------------------------------------------------------------------------------------------------------------------------------------------------------



# Fit the competing model

cohdat$timesince1 <- cohdat$agehl+10
cohdat$timesince2 <- cohdat$agehl+15
cohdat$timesince3 <- cohdat$agehl+25


cohdat$outcomeCompeting <- 1*(cohdat$outcome=="competing")
cohdat$anthra <- 1*(cohdat$anthradosecat!="0 mg/m2")

cohdat <- cohdat[cohdat$exit > cohdat$entry,]

cohdat$id <- 1:nrow(cohdat)
cohdat2 <- tmerge(cohdat,cohdat, id=id, compEvent=event(exit, outcomeCompeting),  tstart=entry,tstop=exit)
cohdat2 <- tmerge(cohdat2,cohdat, id=id, timesince=cumtdc(timesince1, init=0))
cohdat2 <- tmerge(cohdat2,cohdat, id=id, timesince=cumtdc(timesince2))
cohdat2 <- tmerge(cohdat2,cohdat, id=id, timesince=cumtdc(timesince3))


cohdat2$timesince <- factor(cohdat2$timesince, labels=c("5-10yr","10-15yr","15-25yr","25+yr"))


fitCompeting <- coxph(Surv(tstart, tstop, compEvent)~timesince+as.factor(std)+recentsmokedia+ns(txyear, df=3)+relapse5+cardiotoxicRT+I(anthradosecat%in%c("2","3","4")), data=cohdat2, ties="breslow")


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Compute absolute risks for every patient in the case-control data


newdat <- ccdat
subdata <- newdat[,c(1:14,15,19,32:35,39)]
subdata$agehlcat <- cut(subdata$agehl, breaks=c(7,20,30,45), right=F)
subdata$agehlcatnum <- as.numeric(subdata$agehlcat)-1


currentdat <- cbind(subdata[,c("setnr","cc","tumorQuadrant")], subdata[,5:14])
currentdat <- cbind(currentdat, with(subdata, model.matrix(~fam_bc_tot+Menopausal_entry+Menopause_cat_entrynum+AGEFLB_entry+Parity_entry+agehlcatnum+txyearcat))[,-1])

# Compute RR with menopause & parity values at entry

if(mymod == "L1"){
  # L1
  newdat$rr_AR <- c(exp(fitL1$MLE$coef[7:13]%*%t(as.matrix(currentdat[,14:20]))))*
    rowSums((1+exp(fitL1$MLE$coef[2])*currentdat[,4:13])*
              matrix(rep(rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL1$MLE$coef[3:6])), each=2)), nrow(currentdat)), nrow=nrow(currentdat), byrow=TRUE))#/
  
} else if (mymod=="L2"){
  # L2
  newdat$rr_AR <- c(exp(fitL2$MLE$coef[7:13]%*%t(as.matrix(currentdat[,14:20]))))*
    rowSums(exp(fitL2$MLE$coef[2]*currentdat[,4:13])*
              matrix(rep(rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)/sum(rep(exp(c(0,fitL2$MLE$coef[3:6])), each=2)), nrow(currentdat)), nrow=nrow(currentdat), byrow=TRUE))
  
} else if (mymod=="L3"){
  # L3
  newdat$rr_AR <- c(exp(fitL3$MLE$coef[8:14]%*%t(as.matrix(currentdat[,14:20]))))*
    rowSums(sapply(currentdat[,4:13], function(x) exp(model.matrix(~cut(x, breaks=c(0,10,25,500),right=FALSE)-1)[,-1]%*%fitL3$MLE$coef[2:3]))*
              matrix(rep(rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)/sum(rep(exp(c(0,fitL3$MLE$coef[4:7])), each=2)), nrow(currentdat)), nrow=nrow(currentdat), byrow=TRUE))#/
  
  
} else if (mymod=="P1"){
  # P1
  currentdat[,4:13] <- newdat$DMEAN
  newdat$rr_AR <- c(exp(fitP1$MLE$coef[3:9]%*%t(as.matrix(currentdat[,14:20]))))*
    (1+exp(fitP1$MLE$coef[2])*rowMeans(currentdat[,4:13]))
  
} else if (mymod=="P2"){
  # P2
  newdat$agehlcat <- cut(newdat$agehl, breaks=c(7,20,30,45), right=F)
  newdat$agehlcatnum <- as.numeric(newdat$agehlcat)-1
  currentdat <- cbind(newdat[,c("setnr","cc","tumorQuadrant","DMEAN")])
  currentdat <- cbind(currentdat, with(newdat, model.matrix(~fam_bc_tot+Menopausal_entry+Menopause_cat_entrynum+AGEFLB_entry+Parity_entry+agehlcatnum+txyearcat))[,-1])
  newdat$rr_AR <- exp(c(fitP2$MLE$coef[-1]%*%t(as.matrix(currentdat[,4:11]))))
  
} else if (mymod=="P3"){
  # P3
  newdat$agehlcat <- cut(newdat$agehl, breaks=c(7,20,30,45), right=F)
  newdat$agehlcatnum <- as.numeric(newdat$agehlcat)-1
  newdat$DMEANcat <- cut(newdat$DMEAN, breaks=c(0,10,25,500),right=FALSE)
  currentdat <- cbind(newdat[,c("setnr","cc","tumorQuadrant")])
  currentdat <- cbind(currentdat, with(newdat, model.matrix(~DMEANcat+fam_bc_tot+Menopausal_entry+Menopause_cat_entrynum+AGEFLB_entry+Parity_entry+agehlcatnum+txyearcat))[,-1])
  newdat$rr_AR <-exp(c(fitP3$MLE$coef[-1]%*%t(as.matrix(currentdat[,4:12]))))
  
}




newdat$exit <- 100
newdat$compEvent <- 0
newdat$timesince1 <- newdat$agehl+10
newdat$timesince2 <- newdat$agehl+15
newdat$timesince3 <- newdat$agehl+25




newdat$id <- 1:nrow(newdat)
#newdat2 <- tmerge(newdat,newdat, id=id, compEvent=event(exit, outcomeCompeting),  )
newdat2 <- tmerge(newdat,newdat, id=id, timesince=cumtdc(timesince1, init=0), tstart=entry,tstop=exit)
newdat2 <- tmerge(newdat2,newdat, id=id, timesince=cumtdc(timesince2))
newdat2 <- tmerge(newdat2,newdat, id=id, timesince=cumtdc(timesince3))
newdat2$timesince <- factor(newdat2$timesince, labels=c("5-10yr","10-15yr","15-25yr","25+yr")[1:length(unique(newdat2$timesince))])



# Compute all predictions
preds<- lapply(unique(newdat2$id), function(myid){
  
  startage <- newdat[newdat$id==myid,]$entry
  
  stopage1 <- min(newdat[newdat$id==myid,]$entry+1,newdat[newdat$id==myid,]$censage, 75)
  stopage5 <- min(newdat[newdat$id==myid,]$entry+5,newdat[newdat$id==myid,]$censage, 75)
  stopage10 <- min(newdat[newdat$id==myid,]$entry+10,newdat[newdat$id==myid,]$censage, 75)
  stopage20 <- min(newdat[newdat$id==myid,]$entry+20,newdat[newdat$id==myid,]$censage, 75)
  stopage25 <- min(newdat[newdat$id==myid,]$entry+25,newdat[newdat$id==myid,]$censage, 75)
  
  predfit <- survfit(fitCompeting, newdat2[newdat2$id==myid,], id=id, stype=1)
  predfit$time <- predfit$time+newdat[newdat$id==myid,]$entry
  predfit$comphazard <- diff(c(0, predfit$cumhaz))
  
  
  
  
  # Predictions using cohort + general population baseline
  
  tmp <- data.frame(time=predfit$time, comphazard=predfit$comphazard)
  tmp2 <- data.frame(time=liubase$age, bchaz=liubase$baseline*newdat[newdat$id==myid,]$rr_AR)
  
  
  
  
  hazobj2 <- merge(tmp, tmp2, by="time",all=TRUE)
  hazobj2 <- hazobj2[order(hazobj2$time),]
  if(any(is.na(hazobj2$bchaz))) hazobj2[is.na(hazobj2$bchaz),]$bchaz <- 0
  hazobj2[is.na(hazobj2$comphazard),]$comphazard <- 0
  hazobj2 <- hazobj2[hazobj2$time >= startage ,]
  
  hazobj2$Lambdab <- c(0,cumsum(hazobj2$bchaz[-1]))
  hazobj2$Lambdac <- c(0,cumsum(hazobj2$comphazard[-1]))
  
  hazobj2$surv <- exp(-1*hazobj2$Lambdab-hazobj2$Lambdac)
  
  hazobj2$cuminc <- c(0,cumsum((hazobj2$surv*hazobj2$bchaz)[-1]))
  
  
  cuminc_liu_1yr <- sum((hazobj2[hazobj2$time <= stopage1,]$surv*hazobj2[hazobj2$time <= stopage1,]$bchaz)[-1])
  cuminc_liu_5yr <- sum((hazobj2[hazobj2$time <= stopage5,]$surv*hazobj2[hazobj2$time <= stopage5,]$bchaz)[-1])
  cuminc_liu_10yr <- sum((hazobj2[hazobj2$time <= stopage10,]$surv*hazobj2[hazobj2$time <= stopage10,]$bchaz)[-1])
  cuminc_liu_20yr <- sum((hazobj2[hazobj2$time <= stopage20,]$surv*hazobj2[hazobj2$time <= stopage20,]$bchaz)[-1])
  cuminc_liu_25yr <- sum((hazobj2[hazobj2$time <= stopage25,]$surv*hazobj2[hazobj2$time <= stopage25,]$bchaz)[-1])
  
  
  
  
  # Predictions using cohort-based (Langholz-Borgan) baseline
  
  
  tmp3<- data.frame(time=predfit$time, comphazard=predfit$comphazard)
  tmp4 <- data.frame(time=lbbase$age, h1_base=lbbase$h1_base)
  
  hazobj <- merge(tmp3, tmp4, by="time",all=TRUE)
  
  if(any(is.na(hazobj$h1_base))) hazobj[is.na(hazobj$h1_base),]$h1_base <- 0
  if(any(is.na(hazobj$comphazard))) hazobj[is.na(hazobj$comphazard),]$comphazard <- 0
  
  
  
  subhazobj <- hazobj[hazobj$time > startage,]
  
  
  subhazobj$surv <- c(cumprod((1-subhazobj$h1_base*newdat[newdat$id==myid,]$rr_AR-subhazobj$comphazard)))
  subhazobj$cuminc <- cumsum((subhazobj$surv*subhazobj$h1_base*newdat[newdat$id==myid,]$rr_AR))
  
  
  cuminc_lb_1yr <- sum((subhazobj[subhazobj$time <= stopage1,]$h1_base*newdat[newdat$id==myid,]$rr_AR*subhazobj[subhazobj$time <= stopage1,]$surv))
  
  cuminc_lb_5yr <- sum((subhazobj[subhazobj$time <= stopage5,]$h1_base*newdat[newdat$id==myid,]$rr_AR*subhazobj[subhazobj$time <= stopage5,]$surv))
  
  cuminc_lb_10yr <- sum(subhazobj[subhazobj$time <= stopage10,]$h1_base*newdat[newdat$id==myid,]$rr_AR*subhazobj[subhazobj$time <= stopage10,]$surv)
  
  cuminc_lb_20yr <- sum(subhazobj[subhazobj$time <= stopage20,]$h1_base*newdat[newdat$id==myid,]$rr_AR*subhazobj[subhazobj$time <= stopage20,]$surv)
  
  cuminc_lb_25yr <- sum(subhazobj[subhazobj$time <= stopage25,]$h1_base*newdat[newdat$id==myid,]$rr_AR*subhazobj[subhazobj$time <= stopage25,]$surv)
  
  
 
  # Predictions using cohort + survivors in general population baseline
  
  
  
  tmpnkrrho <- data.frame(time=predfit$time, comphazard=predfit$comphazard)
  tmp2nkrrho <- data.frame(time=survbase$age, bchaz=survbase$baseline*newdat[newdat$id==myid,]$rr_AR)
  
  
  
  
  hazobj2nkrrho <- merge(tmpnkrrho, tmp2nkrrho, by="time",all=TRUE)
  hazobj2nkrrho <- hazobj2nkrrho[order(hazobj2nkrrho$time),]
  if(any(is.na(hazobj2nkrrho$bchaz))) hazobj2nkrrho[is.na(hazobj2nkrrho$bchaz),]$bchaz <- 0
  if(any(is.na(hazobj2nkrrho$comphazard))) hazobj2nkrrho[is.na(hazobj2nkrrho$comphazard),]$comphazard <- 0
  hazobj2nkrrho <- hazobj2nkrrho[hazobj2nkrrho$time >= startage ,]
  
  hazobj2nkrrho$Lambdab <- c(0,cumsum(hazobj2nkrrho$bchaz[-1]))
  hazobj2nkrrho$Lambdac <- c(0,cumsum(hazobj2nkrrho$comphazard[-1]))
  
  hazobj2nkrrho$surv <- exp(-1*hazobj2nkrrho$Lambdab-hazobj2nkrrho$Lambdac)
  
  hazobj2nkrrho$cuminc <- c(0,cumsum((hazobj2nkrrho$surv*hazobj2nkrrho$bchaz)[-1]))
  
  
  cuminc_surv_1yr <- sum((hazobj2nkrrho[hazobj2nkrrho$time <= stopage1,]$surv*hazobj2nkrrho[hazobj2nkrrho$time <= stopage1,]$bchaz)[-1])
  cuminc_surv_5yr <- sum((hazobj2nkrrho[hazobj2nkrrho$time <= stopage5,]$surv*hazobj2nkrrho[hazobj2nkrrho$time <= stopage5,]$bchaz)[-1])
  cuminc_surv_10yr <- sum((hazobj2nkrrho[hazobj2nkrrho$time <= stopage10,]$surv*hazobj2nkrrho[hazobj2nkrrho$time <= stopage10,]$bchaz)[-1])
  cuminc_surv_20yr <- sum((hazobj2nkrrho[hazobj2nkrrho$time <= stopage20,]$surv*hazobj2nkrrho[hazobj2nkrrho$time <= stopage20,]$bchaz)[-1])
  cuminc_surv_25yr <- sum((hazobj2nkrrho[hazobj2nkrrho$time <= stopage25,]$surv*hazobj2nkrrho[hazobj2nkrrho$time <= stopage25,]$bchaz)[-1])        
  
  
  data.frame(id=myid,
             cuminc_liu_1yr=cuminc_liu_1yr,
             cuminc_liu_5yr=cuminc_liu_5yr,
             cuminc_liu_10yr=cuminc_liu_10yr,
             cuminc_liu_20yr=cuminc_liu_20yr,
             cuminc_liu_25yr=cuminc_liu_25yr,
             cuminc_lb_1yr=cuminc_lb_1yr,
             cuminc_lb_5yr=cuminc_lb_5yr,
             cuminc_lb_10yr=cuminc_lb_10yr,
             cuminc_lb_20yr=cuminc_lb_20yr,
             cuminc_lb_25yr=cuminc_lb_25yr,
             cuminc_surv_1yr=cuminc_surv_1yr,
             cuminc_surv_5yr=cuminc_surv_5yr,
             cuminc_surv_10yr=cuminc_surv_10yr,
             cuminc_surv_20yr=cuminc_surv_20yr,
             cuminc_surv_25yr=cuminc_surv_25yr
  )
  
})
preds <- as.data.frame(do.call("rbind",preds))



ccdat <- cbind(ccdat, preds)


save(ccdat, file="exampledata_L1preds.RData")



