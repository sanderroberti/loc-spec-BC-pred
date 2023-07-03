library("MASS")
library("survival")
library(Epi)
library(snow)
library(bindata)


load("20221206_rhosim_Smaxcens.RData")

ncores <- 25 # Number of cores to use for parallel computing

simul.data = function(n, m = 2, ifbin, 
                      max.enter = 0, max.t = 50, Smax.t = 0.98, Smax.censor = .90,
                      probz = .5, rxz = 0.1, rxuv = 0.2, rxuz = 0.1, rzuv = 0.2, rzuz = 0.5, 
                      alpha = 1, relx = 2, relz = 1.2, misspecify = F, rho=1, admincens=TRUE){
  
  #### covariates
  Cor = matrix(c(1, rxz, rxuv, rxuz,
                 rxz, 1, rzuv, rzuz,
                 rxuv, rzuv, 1, 0,
                 rxuz, rzuz, 0, 1), ncol = 4)
  
  if( ifbin == T){ Xcov = data.frame(rmvbin(n, c(.5, probz, .5, probz), bincorr = Cor)) }
  if( ifbin == F){ Xcov = data.frame(mvrnorm(n = n, mu = rep(0, 4), Sigma = Cor))
  if(misspecify == T){
    tmpZ = 0.2*Xcov[,1] + Xcov[,4]^4 + rnorm(n,0,1)
    Xcov[,2] = (tmpZ-mean(tmpZ))/sd(tmpZ)
  }
  }
  colnames(Xcov) = c("X", "Z", "Uv", "Uz")
  
  #### coefficients
  beta0 = log(rho*-log(Smax.t))
  beta1 = log(relx)
  beta2 = log(relz)
  lam = with(Xcov, exp((beta0 + beta1*X + beta2*Z)/alpha))
  
  #### survival data
  survtime = mapply(rweibull, n = 1, shape = alpha, scale = 1/lam)
  entrtime = runif(n, 0, max.enter)
  deadtime = rexp(n, rate = -log(Smax.censor))
  if(admincens){
    censtime <- max.t - entrtime
  } else{
    censtime = 9999999
  }
  
  eventime = mapply(min, survtime, censtime, deadtime)
  ind.fail = (survtime <= mapply(min, censtime, deadtime))
  outcome <- ifelse(ind.fail,"event", ifelse(censtime<deadtime, "admin","comp"))
  
  DATA = data.frame(cbind(Xcov, eventime, ind.fail, outcome))
  
  
  
  #### MCAR with n(MAR/NCC) ####
  dat_mcar = DATA
  
  
  return(list(mcar = dat_mcar))#, mar_cc = dat_mar_cc, mar_ncc = dat_mar_ncc))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(82849583)


results <- NULL

for(rrX in c(1,2,4)){
  for(rrZ in c(1,2,4)){
    for(myrho in c(1,2,7,10)){ # 
      
      mysmaxt <- .98
      mymaxt <- 15
      mysmaxcens <- exp(Smaxres[Smaxres$rrX==rrX & Smaxres$rrZ==rrZ & Smaxres$rho==myrho,]$logSmaxcens2)
      
      coh_base <- simul.data(n = 2000000, ifbin = T, relx = rrX, relz = rrZ, misspecify = F,rxz=.6, rxuv=0, rxuz=0, rzuv=0, rzuz=0, Smax.t=mysmaxt, Smax.censor=mysmaxcens, max.t=mymaxt)$mcar
      coh_base$dum <- 1
      coh_base$entry <- 0
      
      cph <- coxph(Surv(entry,eventime, ind.fail)~dum, data=coh_base)
      cph$coef["dum"] <- 0

      for(N in c(10000, 20000, 30000, 50000, 100000)){ # 
        print(paste(rrX, rrZ, myrho, N))
        cl <- makeCluster(ncores, outfile="rhosim_part1_debug.txt")
        clusterEvalQ(cl=cl, {RNGkind("L'Ecuyer-CMRG"); library(survival);library(Epi);library(MASS); library(RhpcBLASctl);library(bindata); blas_set_num_threads(1)})
        clusterExport(cl, c("simul.data","N","rrX","rrZ","coh_base","myrho","cph","mysmaxcens","mysmaxt","mymaxt"))
        
        result_list <- clusterMap(cl=cl,1:500, fun=function(it){
          
          # Generate separate cohorts for modeling and to act as registry
          coh_model <- simul.data(n = N, ifbin = T, relx = rrX, relz = rrZ, misspecify = F, rho=myrho, Smax.t=mysmaxt, rxz=.6, rxuv=0, rxuz=0, rzuv=0, rzuz=0, Smax.censor=mysmaxcens, max.t=mymaxt)$mcar #.45
          coh_model$dum <- 1
          
          censrate <- 1-mean(coh_model$ind.fail)
          admincens <- mean(coh_model$outcome=="admin")
          
         
          coh_model$entry <- 0
  
          
          coh_modelclone1 <- coh_model
          coh_modelclone1$eventime <- pmin(median(coh_model[coh_model$ind.fail==T,]$eventime), coh_modelclone1$eventime)
          mypreds1 <- predict(cph,newdata=coh_modelclone1, type="expected")
          mypreds1[is.na(mypreds1)] <- 0

          coh_modelclone2 <- coh_model
          coh_modelclone2$entry <- median(coh_model[coh_model$ind.fail==T,]$eventime)
          mypreds2 <- predict(cph,newdata=coh_modelclone2, type="expected")
          mypreds2[is.na(mypreds2)] <- 0
          
          coh_model$popinc1 <- mypreds1
          coh_model$popinc2 <- mypreds2
          
          m <- 3
          
          coh_model$nrisk = with(coh_model, mapply(function(i){ sum(eventime >= eventime[i]) }, 1:nrow(coh_model)))
          coh_model$incl.prob = with(coh_model, mapply(function(i){
            1-prod(1-m/(nrisk[which((eventime < eventime[i]) & (ind.fail == 1))]-1))
          }, 1:nrow(coh_model)))
          coh_model$samuel.wgt = 1/coh_model$incl.prob
          
          # Function to compute incidence by individual ages
          incfun <- function(ages){
            sapply(ages, function(age){
              sum(coh_base[floor(coh_base$eventime)==age,]$ind.fail)/
                sum(floor(coh_base$eventime)>=age)
            })
          }
          
          # Case-control sampling using Epi package
          ncc <- ccwc(exit=eventime,fail=ind.fail, data=coh_model,include=list(X,Z, eventime, samuel.wgt, nrisk, popinc1, popinc2), controls=m,  silent=TRUE)
          
          fit.ncc <- clogit(Fail~X+Z+strata(Set), data=ncc)
          
          rr_AR <- exp(as.matrix(ncc[,c("X","Z")])%*%fit.ncc$coefficients)
          
          
          fit.coh <- coxph(Surv(eventime,ind.fail) ~X+Z, data=coh_model)
          coh_model$RR <- exp(predict(fit.coh, type="lp"))
          
          coh_basehaz <- basehaz(fit.coh, centered=FALSE)
          
          ARcoh <- 1-1/mean(coh_model$RR)
          AR1coh <- 1-1/mean(coh_model[coh_model$eventime<median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          AR2coh <- 1-1/mean(coh_model[coh_model$eventime>=median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          
          ARcohbruzzi <- 1-mean(1/coh_model[coh_model$ind.fail==TRUE,]$RR)
          AR1cohbruzzi <- 1-mean(1/coh_model[coh_model$ind.fail==TRUE & coh_model$eventime<median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          AR2cohbruzzi <- 1-mean(1/coh_model[coh_model$ind.fail==TRUE & coh_model$eventime>=median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          
          
          
          # Rho computation based on case-control study
          res2 <- ncc[,c("eventime","Set","Fail","samuel.wgt","popinc1","popinc2")]
          res2 <- cbind(res2, RR=rr_AR)
          res2 <- res2[order(res2$eventime),]
          res2$natrisk <- sapply(res2$eventime, function(x) sum(coh_model$eventime >= x))
          res2$setsize <- sapply(res2$Set, function(xset) sum(ncc$Set==xset))
          res2$wgt_samuelsen <- NA
          for (i in 1:nrow(res2)){
            if (res2$Fail[i] == 1){
              res2$wgt_samuelsen[i] <- 1
            } else {
              res2$wgt_samuelsen[i] <- 1/(1-prod(1-(res2[res2$eventime < res2$eventime[i] & res2$Fail==1,]$setsize-1)/(res2[res2$eventime < res2$eventime[i] & res2$Fail==1,]$natrisk-1)))
            }
          }
          
          res2[res2$Fail==1,]$samuel.wgt <- 1
          
          res2[res2$Fail==0,]$wgt_samuelsen <- sum(1-coh_model$ind.fail)*res2[res2$Fail==0,]$wgt_samuelsen/sum(res2[res2$Fail==0,]$wgt_samuelsen)
          
          ARncc <- 1-1/(sum(res2$RR*res2$wgt_samuelsen)/sum(res2$wgt_samuelsen))
          AR1ncc <- 1-1/(sum(res2[res2$eventime < median(ncc[ncc$Fail==1,]$eventime),]$RR*res2[res2$eventime < median(ncc[ncc$Fail==1,]$eventime),]$wgt_samuelsen)/sum(res2[res2$eventime < median(ncc[ncc$Fail==1,]$eventime),]$wgt_samuelsen))
          AR2ncc <- 1-1/(sum(res2[res2$eventime >= median(ncc[ncc$Fail==1,]$eventime),]$RR*res2[res2$eventime >= median(ncc[ncc$Fail==1,]$eventime),]$wgt_samuelsen)/sum(res2[res2$eventime >= median(ncc[ncc$Fail==1,]$eventime),]$wgt_samuelsen))
          
          ARnccbruzzi <- 1-mean(1/res2[res2$Fail==1,]$RR)
          AR1nccbruzzi <- 1-mean(1/res2[res2$Fail==1 & res2$eventime < median(ncc[ncc$Fail==1,]$eventime),]$RR)
          AR2nccbruzzi <- 1-mean(1/res2[res2$Fail==1 & res2$eventime >= median(ncc[ncc$Fail==1,]$eventime),]$RR)
          
          
          # rho_ncc is computed using time-dependent phi (cut at median event time)
          # rho_ncc0 is computed using overall phi
          # rho_nccbruzzi and rho_ncc0bruzzi are computed using time-dependent and overall phi from Bruzzi formula, respectively
          
          rho_ncc0 <- sum(ncc$Fail)/sum(res2$RR*(1-ARncc)*(res2$popinc1+res2$popinc2)*res2$wgt_samuelsen)
          rho_ncc <- sum(ncc$Fail)/sum(res2$RR*((1-AR1ncc)*res2$popinc1+(1-AR2ncc)*res2$popinc2)*res2$wgt_samuelsen)
          
          rho_ncc0bruzzi <- sum(ncc$Fail)/sum(res2$RR*(1-ARnccbruzzi)*(res2$popinc1+res2$popinc2)*res2$wgt_samuelsen)
          rho_nccbruzzi <- sum(ncc$Fail)/sum(res2$RR*((1-AR1nccbruzzi)*res2$popinc1+(1-AR2nccbruzzi)*res2$popinc2)*res2$wgt_samuelsen)
          
          
          # Rho based on cohort
          # rho_coh is computed using time-dependent phi (cut at median event time)
          # rho_coh0 is computed using overall phi
          # rho_cohbruzzi and rho_coh0bruzzi are computed using time-dependent and overall phi from Bruzzi formula, respectively
         
          
          rho_coh <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1*(1-AR1coh)+coh_model$popinc2*(1-AR2coh)))
          rho_coh0 <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1+coh_model$popinc2)*(1-ARcoh))
          
          rho_cohbruzzi <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1*(1-AR1cohbruzzi)+coh_model$popinc2*(1-AR2cohbruzzi)))
          rho_coh0bruzzi <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1+coh_model$popinc2)*(1-ARcohbruzzi))
          
          data.frame(N=N, it=it,rrX=rrX,rrZ=rrZ,rho=myrho,censrate=censrate,admincens=admincens, rho_ncc=rho_ncc,rho_ncc0=rho_ncc0, rho_coh=rho_coh,rho_coh0=rho_coh0,rho_nccbruzzi=rho_nccbruzzi,rho_ncc0bruzzi=rho_ncc0bruzzi, rho_cohbruzzi=rho_cohbruzzi,rho_coh0bruzzi=rho_coh0bruzzi, ARncc=ARncc,AR1ncc=AR1ncc,AR2ncc=AR2ncc,ARcoh=ARcoh,AR1coh=AR1coh,AR2coh=AR2coh, ARnccbruzzi=ARnccbruzzi,AR1nccbruzzi=AR1nccbruzzi,AR2nccbruzzi=AR2nccbruzzi,ARcohbruzzi=ARcohbruzzi,AR1cohbruzzi=AR1cohbruzzi,AR2cohbruzzi=AR2cohbruzzi, logrrX_ncc=fit.ncc$coefficients[1], logrrZ_ncc=fit.ncc$coefficients[2], logrrX_coh=fit.coh$coefficients[1], logrrZ_coh=fit.coh$coefficients[2])
          
        })
        
        
        tmp <- do.call(rbind, result_list)
        results <- rbind(results, tmp)
        closeAllConnections()
      }
    }
  }
}


save(results, file="rhosim_part1_results.RData")
