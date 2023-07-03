library("MASS")
library("survival")
library(Epi)
library(snow)
library(bindata)
library(survey)

load("20230102_ageeffsim_Smaxcens.RData")

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
  
  Xcov$X <- 1*(Xcov$X>.84)
  
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
  
  
  return(list(mcar = dat_mcar))
}




RNGkind("L'Ecuyer-CMRG")
set.seed(82849511)


results <- NULL

for(rrX in c(.5, 1, 1.5)){
  for(rrZ in c(1,1.25,1.5,2)){
    for(myrho in c(1, 2, 7, 10)){
      mysmaxt <- .98
      mymaxt <- 15
      mysmaxcens <- exp(Smaxres[Smaxres$rrX==rrX & Smaxres$rrZ==rrZ & Smaxres$rho==myrho,]$logSmaxcens2)
      
      
      coh_base <- simul.data(n = 2000000, ifbin = F, relx = rrX, relz = rrZ, rxz=0.6, misspecify = F, Smax.t=mysmaxt, rxuv=0, rxuz=0, rzuv=0, rzuz=0, Smax.censor=mysmaxcens, max.t=mymaxt)$mcar
      
      censratebasecoh <- 1-mean(coh_base$ind.fail)
      coh_base$Zcat2 <- cut(coh_base$Z, breaks=qnorm(0:10/10))
      
      coh_base$dum <- 1
      coh_base$entry <- 0
      
      cph <- coxph(Surv(entry,eventime, ind.fail)~dum, data=coh_base, control = coxph.control(timefix = FALSE))
      cph$coef["dum"] <- 0
      
      
      for(N in c(10000, 50000, 100000)){
        print(paste(rrX, rrZ, myrho, N))
        
        cl <- makeCluster(25, outfile="ageeffsim_debug.txt")
        clusterEvalQ(cl=cl, {RNGkind("L'Ecuyer-CMRG"); library(survival);library(Epi);library(MASS); library(RhpcBLASctl);library(bindata);library(survey); blas_set_num_threads(1)})
        clusterExport(cl, c("simul.data","N","rrX","rrZ","myrho","coh_base","cph","mysmaxt","mysmaxcens","mymaxt"))
        
        result_list <- clusterMap(cl=cl,1:500, fun=function(it){
          
          
          
          
          coh_model <- simul.data(n = N, ifbin = F, relx = rrX, relz = rrZ, rxz=0.6, misspecify = F, rho=myrho, Smax.t=mysmaxt, rxuv=0, rxuz=0, rzuv=0, rzuz=0, Smax.censor=mysmaxcens, max.t=mymaxt)$mcar
          
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
          
          coh_model$Zcat <- cut(coh_model$Z, breaks=qnorm(0:10/10))
          coh_model$Zcat2 <- cut(coh_model$Z, breaks=qnorm(0:10/10))
          
          m <- 3
          
          coh_model$nrisk = with(coh_model, mapply(function(i){ sum(eventime >= eventime[i]) }, 1:nrow(coh_model)))
          coh_model$incl.prob = with(coh_model, mapply(function(i){
            1-prod(1-m/(nrisk[which((eventime < eventime[i]) & (ind.fail == 1))]-1))
          }, 1:nrow(coh_model)))
          coh_model$samuel.wgt = 1/coh_model$incl.prob
          
          
          
          
          coh_model$nriskstrat <- NA
          coh_model$incl.probstrat <- NA
          coh_model$samuel.wgtstrat <- NA
          
          for(mystrat in unique(coh_model$Zcat2)){
            
            coh_model[coh_model$Zcat2==mystrat,]$nriskstrat <- with(coh_model[coh_model$Zcat2==mystrat,], mapply(function(i){ sum(eventime >= eventime[i]) }, 1:nrow(coh_model[coh_model$Zcat2==mystrat,])))
            
          }
          
          
          
          
          # Case-control sampling using Epi package
          ncc <- ccwc(exit=eventime,fail=ind.fail, data=coh_model,include=list(X,Z,Zcat2, eventime, nriskstrat, popinc1, popinc2), controls=m, match=list(Zcat), silent=TRUE)
          
          ncc$mstrat <- sapply(1:nrow(ncc), function(myrow) sum(ncc$Zcat2==ncc[myrow,]$Zcat2 & ncc$Set==ncc[myrow,]$Set & ncc$Fail==0))
          
          ncc$incl.probstrat <- with(ncc, mapply(function(i){
            1-prod(1-mstrat[i]/(nriskstrat[which((eventime < eventime[i]) & (Fail == 1) & Zcat2==Zcat2[i])]-1))
          }, 1:nrow(ncc)))
          
          ncc$samuel.wgtstrat = 1/ncc$incl.probstrat
          
          ncc[ncc$Fail==1,]$samuel.wgtstrat<- 1
          
          ncc$samwgt2 <- ncc$samuel.wgtstrat
          
          for(mystrat in unique(ncc$Zcat2)){
            ncc[ncc$Zcat2==mystrat & ncc$Fail==0,]$samwgt2 <- ncc[ncc$Zcat2==mystrat & ncc$Fail==0,]$samwgt2*sum(coh_model$Zcat2==mystrat & coh_model$ind.fail==FALSE)/sum(ncc[ncc$Zcat2==mystrat & ncc$Fail==0,]$samwgt2)
          }
          
          
          dsgn <- svydesign(ids=~1,strata=~Zcat2, data=ncc, weights=ncc$samwgt2)
          
          ccfit <- tryCatch(svyglm(Fail~X+Z, data=ncc, family="binomial", design=dsgn), error=function(e) list(coefficients=rep(NA,3)))
          
          cohfit <- tryCatch(coxph(Surv(eventime,ind.fail)~X+Z, data=coh_model, control = coxph.control(timefix = FALSE)), error=function(e) list(coefficients=rep(NA,2)))
          
          
          
          
          rr_AR <- exp(as.matrix(ncc[,c("X","Z")])%*%ccfit$coefficients[2:3])
          
          
          coh_model$RR <- exp(predict(cohfit, type="lp"))
          
          
          ARcoh <- 1-1/mean(coh_model$RR)
          AR1coh <- 1-1/mean(coh_model[coh_model$eventime<median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          AR2coh <- 1-1/mean(coh_model[coh_model$eventime>=median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          
          ARcohbruzzi <- 1-mean(1/coh_model[coh_model$ind.fail==TRUE,]$RR)
          AR1cohbruzzi <- 1-mean(1/coh_model[coh_model$ind.fail==TRUE & coh_model$eventime<median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          AR2cohbruzzi <- 1-mean(1/coh_model[coh_model$ind.fail==TRUE & coh_model$eventime>=median(coh_model[coh_model$ind.fail==T,]$eventime),]$RR)
          
          
          
          ARncc <- 1-1/(sum(rr_AR*ncc$samwgt2)/sum(ncc$samwgt2))
          AR1ncc <- 1-1/(sum(rr_AR[ncc$eventime < median(coh_model[coh_model$ind.fail==T,]$eventime)]*ncc[ncc$eventime < median(coh_model[coh_model$ind.fail==T,]$eventime),]$samwgt2)/sum(ncc[ncc$eventime < median(coh_model[coh_model$ind.fail==T,]$eventime),]$samwgt2))
          AR2ncc <- 1-1/(sum(rr_AR[ncc$eventime >= median(coh_model[coh_model$ind.fail==T,]$eventime)]*ncc[ncc$eventime >= median(coh_model[coh_model$ind.fail==T,]$eventime),]$samwgt2)/sum(ncc[ncc$eventime >=median(coh_model[coh_model$ind.fail==T,]$eventime),]$samwgt2))
          
          ARnccbruzzi <- 1-mean(1/rr_AR[ncc$Fail==1])
          AR1nccbruzzi <- 1-mean(1/rr_AR[ncc$Fail==1 & ncc$eventime < median(coh_model[coh_model$ind.fail==T,]$eventime)])
          AR2nccbruzzi <- 1-mean(1/rr_AR[ncc$Fail==1 & ncc$eventime >= median(coh_model[coh_model$ind.fail==T,]$eventime)])
          
          
          # rho_ncc is computed using time-dependent phi (cut at median event time)
          # rho_ncc0 is computed using overall phi
          # rho_nccbruzzi and rho_ncc0bruzzi are computed using time-dependent and overall phi from Bruzzi formula, respectively
          
          rho_ncc0 <- sum(ncc$Fail)/sum(rr_AR*(1-ARncc)*(ncc$popinc1+ncc$popinc2)*ncc$samwgt2)
          rho_ncc <- sum(ncc$Fail)/sum(rr_AR*((1-AR1ncc)*ncc$popinc1+(1-AR2ncc)*ncc$popinc2)*ncc$samwgt2)
          
          rho_ncc0bruzzi <- sum(ncc$Fail)/sum(rr_AR*(1-ARnccbruzzi)*(ncc$popinc1+ncc$popinc2)*ncc$samwgt2)
          rho_nccbruzzi <- sum(ncc$Fail)/sum(rr_AR*((1-AR1nccbruzzi)*ncc$popinc1+(1-AR2nccbruzzi)*ncc$popinc2)*ncc$samwgt2)
          
          
          # Rho based on cohort
          # rho_coh is computed using time-dependent phi (cut at median event time)
          # rho_coh0 is computed using overall phi
          # rho_cohbruzzi and rho_coh0bruzzi are computed using time-dependent and overall phi from Bruzzi formula, respectively
          
          rho_coh <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1*(1-AR1coh)+coh_model$popinc2*(1-AR2coh)))
          rho_coh0 <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1+coh_model$popinc2)*(1-ARcoh))
          
          rho_cohbruzzi <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1*(1-AR1cohbruzzi)+coh_model$popinc2*(1-AR2cohbruzzi)))
          rho_coh0bruzzi <- sum(coh_model$ind.fail)/sum(coh_model$RR*(coh_model$popinc1+coh_model$popinc2)*(1-ARcohbruzzi))
          
          
          
          
          list(it=it, censrate=censrate,admincens=admincens, ccfit=ccfit$coefficients, cohfit=cohfit$coefficients, rho_ncc=rho_ncc,rho_ncc0=rho_ncc0, rho_coh=rho_coh,rho_coh0=rho_coh0,rho_nccbruzzi=rho_nccbruzzi,rho_ncc0bruzzi=rho_ncc0bruzzi, rho_cohbruzzi=rho_cohbruzzi,rho_coh0bruzzi=rho_coh0bruzzi, ARncc=ARncc,AR1ncc=AR1ncc,AR2ncc=AR2ncc,ARcoh=ARcoh,AR1coh=AR1coh,AR2coh=AR2coh, ARnccbruzzi=ARnccbruzzi,AR1nccbruzzi=AR1nccbruzzi,AR2nccbruzzi=AR2nccbruzzi,ARcohbruzzi=ARcohbruzzi,AR1cohbruzzi=AR1cohbruzzi,AR2cohbruzzi=AR2cohbruzzi)
          
        })

        results <- c(results, list(list(N=N,rrX=rrX,rrZ=rrZ,rho=myrho,censratebasecoh=censratebasecoh, results=result_list)))
        
        closeAllConnections()
      }
    }
  }
}

save(results, file="ageeffsim_results_full.RData")
