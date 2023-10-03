library(bbmle)
library(car)
library(numDeriv)
library(data.table)
library(Matrix)




linERRloglikcomposite <- function(params, data, doses, set, status, loc, in_cohort, corrvars=NULL, corrvars_exclude_conditional=NULL, ccmethod="CCAL", weights){ # -log(likelihood)
  
  intc <- params[1]
  beta0 <- params[2]
  
  if(ccmethod %in% c("CCAL","CL")){
    alpha <- c(0,0,rep(params[3:(length(doses)/2+1)],each=2))
  } else {
    alpha <- rep(0,length(doses)/2)
  }
  
  delta <- tail(params, length(corrvars))
  
  wgtdata <- data[data[,in_cohort]==1,]
  
  if(!ccmethod=="meandose"){
    
    
    if(ccmethod=="CCAL"){
      # CCAL
      tmp <- exp(as.matrix(wgtdata[,corrvars,drop=FALSE])%*%delta+intc)
      pmat <- matrix(rep(tmp, length(doses)), ncol=length(doses))
      tmp2 <- (1+beta0*wgtdata[,doses])
      tmp3 <- tmp2*matrix(rep(exp(alpha), nrow(wgtdata)), ncol=length(doses), byrow=TRUE)
      pmat <- pmat*tmp3
      pmat <- cbind(1, pmat)
      pmat <- as.data.table(pmat)
      rsums <- rowSums(pmat)
      pmat <- pmat[, (names(pmat)) := lapply(.SD,function(x) log(x/rsums ))]
      
      
      Ymat <- matrix(0, nrow=nrow(wgtdata), ncol=length(doses)+1)
      Ymat[cbind(1:nrow(wgtdata), wgtdata[,loc]+1)] <- 1
      Ymat[wgtdata[,status]==0,1] <- 1
      Ymat[wgtdata[,status]==0,-1] <- 0
      
      minuslogl1 <- -1*sum(rowSums(pmat*Ymat)*wgtdata[,weights])
      
    } 
  } else {
    
    tmp <- exp(as.matrix(wgtdata[,corrvars,drop=FALSE])%*%delta+intc)
    
    tmp2 <- (1+beta0*rowMeans(wgtdata[,doses]))
    
    pmat <- tmp*tmp2
    pmat <- pmat/(1+pmat)
    
    minuslogl1 <- -1*sum(wgtdata[,weights]*(wgtdata[,status]*log(pmat)+(1-wgtdata[,status])*log(1-pmat)))
    
  }
  
  
  
  
  conddata <- data[data[,set] %in% unique(data[data[,status]==1 & data[,in_cohort]==0, set]),]
  
  conddata <- conddata[conddata[,set] %in% as.numeric(names(which(table(conddata[,set],conddata[,status])[,2]==1))),]
  if(ccmethod!="CL") conddata <- conddata[conddata[,set] %in% as.numeric(names(which(table(conddata[,set])>1))),]
  
  
  delta2 <- delta[-corrvars_exclude_conditional]
  corrvars2 <- corrvars[-corrvars_exclude_conditional]
  
  if(!ccmethod=="meandose"){
    
    rrmat <- matrix(rep(exp(as.matrix(conddata[,corrvars2,drop=FALSE])%*%delta2), length(doses)), ncol=length(doses))*(1+beta0*conddata[,doses])*matrix(rep(exp(alpha), nrow(conddata)), ncol=length(doses), byrow=TRUE)
    
    if(ccmethod=="CCAL"){
      # CCAL
      
      rrcases <- rrmat[conddata[,status]==1,][cbind(1:sum(conddata[,status]), conddata[conddata[,status]==1,loc])]
      
      rrdf <- data.frame(set=conddata[,set], id=1,rrmat)
      rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
      rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")
      
      minuslogl2 <- -(sum(log(rrcases)) - sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    }
  } else {
    rrmat <- (1+beta0*rowMeans(conddata[,doses]))*exp(as.matrix(conddata[,corrvars2,drop=FALSE])%*%delta2)
    rrdf <- data.frame(set=conddata[,set], id=1, rrmat)
    rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
    rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")
    
    minuslogl2 <- -(sum(log(rrmat[conddata[,status]==1,]))-sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    
  }
  
  
  minuslogl <- minuslogl1+minuslogl2
  
  minuslogl
}


linERRloglikcompositeexp <- function(params, data, doses, set, status, loc, in_cohort, corrvars=NULL, corrvars_exclude_conditional=NULL, ccmethod="CCAL", weights){ # -log(likelihood)
  
  intc <- params[1]
  beta0 <- params[2]
  
  if(ccmethod %in% c("CCAL","CL")){
    alpha <- c(0,0,rep(params[3:(length(doses)/2+1)],each=2))
  } else {
    alpha <- rep(0,length(doses)/2)
  }
  
  delta <- tail(params, length(corrvars))
  
  wgtdata <- data[data[,in_cohort]==1,]
  
  if(!ccmethod=="meandose"){
    
    
    if(ccmethod=="CCAL"){
      # CCAL
      tmp <- exp(as.matrix(wgtdata[,corrvars,drop=FALSE])%*%delta+intc)
      pmat <- matrix(rep(tmp, length(doses)), ncol=length(doses))
      tmp2 <- exp(beta0*wgtdata[,doses])
      tmp3 <- tmp2*matrix(rep(exp(alpha), nrow(wgtdata)), ncol=length(doses), byrow=TRUE)
      pmat <- pmat*tmp3
      pmat <- cbind(1, pmat)
      pmat <- as.data.table(pmat)
      rsums <- rowSums(pmat)
      pmat <- pmat[, (names(pmat)) := lapply(.SD,function(x) log(x/rsums ))]
      
      
      Ymat <- matrix(0, nrow=nrow(wgtdata), ncol=length(doses)+1)
      Ymat[cbind(1:nrow(wgtdata), wgtdata[,loc]+1)] <- 1
      Ymat[wgtdata[,status]==0,1] <- 1
      Ymat[wgtdata[,status]==0,-1] <- 0
      
      minuslogl1 <- -1*sum(rowSums(pmat*Ymat)*wgtdata[,weights])
      
    } 
  } else {
    
    tmp <- exp(as.matrix(wgtdata[,corrvars,drop=FALSE])%*%delta+intc)
    
    tmp2 <- exp(beta0*rowMeans(wgtdata[,doses]))
    
    pmat <- tmp*tmp2
    pmat <- pmat/(1+pmat)
    
    minuslogl1 <- -1*sum(wgtdata[,weights]*(wgtdata[,status]*log(pmat)+(1-wgtdata[,status])*log(1-pmat)))
    
  }
  
  
  
  
  conddata <- data[data[,set] %in% unique(data[data[,status]==1 & data[,in_cohort]==0, set]),]
  
  conddata <- conddata[conddata[,set] %in% as.numeric(names(which(table(conddata[,set],conddata[,status])[,2]==1))),]
  if(ccmethod!="CL") conddata <- conddata[conddata[,set] %in% as.numeric(names(which(table(conddata[,set])>1))),]
  
  
  delta2 <- delta[-corrvars_exclude_conditional]
  corrvars2 <- corrvars[-corrvars_exclude_conditional]
  
  if(!ccmethod=="meandose"){
    
    rrmat <- matrix(rep(exp(as.matrix(conddata[,corrvars2,drop=FALSE])%*%delta2), length(doses)), ncol=length(doses))*(exp(beta0*conddata[,doses]))*matrix(rep(exp(alpha), nrow(conddata)), ncol=length(doses), byrow=TRUE)
    
    if(ccmethod=="CCAL"){
      # CCAL
      
      rrcases <- rrmat[conddata[,status]==1,][cbind(1:sum(conddata[,status]), conddata[conddata[,status]==1,loc])]
      
      rrdf <- data.frame(set=conddata[,set], id=1,rrmat)
      rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
      rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")
      
      minuslogl2 <- -(sum(log(rrcases)) - sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    }
  } else {
    rrmat <- exp(beta0*rowMeans(conddata[,doses]))*exp(as.matrix(conddata[,corrvars2,drop=FALSE])%*%delta2)
    rrdf <- data.frame(set=conddata[,set], id=1, rrmat)
    rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
    rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")
    
    minuslogl2 <- -(sum(log(rrmat[conddata[,status]==1,]))-sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    
  }
  
  
  minuslogl <- minuslogl1+minuslogl2
  
  minuslogl
}


linERRloglikcompositecat <- function(params, data, doses, set, status, loc, in_cohort, corrvars=NULL, corrvars_exclude_conditional=NULL, ccmethod="CCAL", weights){ # -log(likelihood)
  
  intc <- params[1]
  beta0 <- params[2:3]
  
  if(ccmethod %in% c("CCAL","CL")){
    alpha <- c(0,0,rep(params[4:(length(doses)/2+2)],each=2))
  } else {
    alpha <- rep(0,length(doses)/2)
  }
  
  delta <- tail(params, length(corrvars))
  
  wgtdata <- data[data[,in_cohort]==1,]
  
  if(!ccmethod=="meandose"){
    
    
    if(ccmethod=="CCAL"){
      # CCAL
      tmp <- exp(as.matrix(wgtdata[,corrvars,drop=FALSE])%*%delta+intc)
      pmat <- matrix(rep(tmp, length(doses)), ncol=length(doses))
      tmp2 <- sapply(wgtdata[,doses], function(x) exp(model.matrix(~cut(x, breaks=c(0,10,25,500),right=FALSE)-1)[,-1]%*%beta0))
      tmp3 <- tmp2*matrix(rep(exp(alpha), nrow(wgtdata)), ncol=length(doses), byrow=TRUE)
      
      pmat <- pmat*tmp3
      pmat <- cbind(1, pmat)
      pmat <- as.data.table(pmat)
      rsums <- rowSums(pmat)
      pmat <- pmat[, (names(pmat)) := lapply(.SD,function(x) log(x/rsums ))]
      
      
      Ymat <- matrix(0, nrow=nrow(wgtdata), ncol=length(doses)+1)
      Ymat[cbind(1:nrow(wgtdata), wgtdata[,loc]+1)] <- 1
      Ymat[wgtdata[,status]==0,1] <- 1
      Ymat[wgtdata[,status]==0,-1] <- 0
      
      minuslogl1 <- -1*sum(rowSums(pmat*Ymat)*wgtdata[,weights])
      
    } 
  } else {
    
    tmp <- exp(as.matrix(wgtdata[,corrvars,drop=FALSE])%*%delta+intc)
    
    tmp2 <- exp(model.matrix(~cut(rowMeans(wgtdata[,doses]), breaks=c(0,10,25,500), right=FALSE)-1)[,-1]%*%beta0)
    
    pmat <- tmp*tmp2
    pmat <- pmat/(1+pmat)
    
    minuslogl1 <- -1*sum(wgtdata[,weights]*(wgtdata[,status]*log(pmat)+(1-wgtdata[,status])*log(1-pmat)))
    
  }
  
  
  
  
  conddata <- data[data[,set] %in% unique(data[data[,status]==1 & data[,in_cohort]==0, set]),]
  
  conddata <- conddata[conddata[,set] %in% as.numeric(names(which(table(conddata[,set],conddata[,status])[,2]==1))),]
  if(ccmethod!="CL") conddata <- conddata[conddata[,set] %in% as.numeric(names(which(table(conddata[,set])>1))),]
  
  
  delta2 <- delta[-corrvars_exclude_conditional]
  corrvars2 <- corrvars[-corrvars_exclude_conditional]
  
  if(!ccmethod=="meandose"){
    
    rrmat <- matrix(rep(exp(as.matrix(conddata[,corrvars2,drop=FALSE])%*%delta2), length(doses)), ncol=length(doses))*sapply(conddata[,doses], function(x) exp(model.matrix(~cut(x, breaks=c(0,10,25,500),right=FALSE)-1)[,-1]%*%beta0))*matrix(rep(exp(alpha), nrow(conddata)), ncol=length(doses), byrow=TRUE)
    
    if(ccmethod=="CCAL"){
      # CCAL
      
      rrcases <- rrmat[conddata[,status]==1,][cbind(1:sum(conddata[,status]), conddata[conddata[,status]==1,loc])]
      
      rrdf <- data.frame(set=conddata[,set], id=1,rrmat)
      rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
      rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")
      
      minuslogl2 <- -(sum(log(rrcases)) - sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    }
  } else {
    rrmat <- exp(model.matrix(~cut(rowMeans(conddata[,doses]), breaks=c(0,10,25,500), right=FALSE)-1)[,-1]%*%beta0)*exp(as.matrix(conddata[,corrvars2,drop=FALSE])%*%delta2)
    rrdf <- data.frame(set=conddata[,set], id=1, rrmat)
    rrdf$id <- with(rrdf, ave(as.numeric(set), set, FUN = seq_along))
    rrdfwide <- reshape(rrdf, direction="wide", idvar="set",timevar="id")
    
    minuslogl2 <- -(sum(log(rrmat[conddata[,status]==1,]))-sum(log(rowSums(rrdfwide[,-1], na.rm=TRUE))))
    
  }
  
  
  minuslogl <- minuslogl1+minuslogl2
  
  minuslogl
}


linearERRfitcomposite <- function(data, doses, set, status, loc, in_cohort, corrvars=NULL,corrvars_exclude_conditional=NULL,repar=FALSE, ccmethod="CCAL", initpars=rep(0,1+length(doses)/2+length(corrvars)), fitopt=list(maxit=5000), fitNull=TRUE, useOld=FALSE, uplimBeta=5, weights){
  
  if(ccmethod=="CL") corrvars <- NULL
  
  opt_method <- ifelse(repar, "Nelder-Mead", "L-BFGS-B")
  
  
  if(is.null(fitopt$maxit)){
    fitopt <- c(fitopt, list(maxit=5000))
  }
  if(is.null(fitopt$reltol) & opt_method!="L-BFGS-B"){
    fitopt <- c(fitopt, list(reltol=1e-10))
  }
  if (is.null(fitopt$pgtol) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(pgtol=1e-8))
  }
  if (is.null(fitopt$factr) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(factr=1e4))
  }
  if (is.null(fitopt$ndeps) & opt_method=="L-BFGS-B"){
    parlen <- 2+ifelse(ccmethod %in% c("CCAL","CL"),length(doses)/2-1,0)+length(corrvars)
    fitopt <- c(fitopt, list(ndeps=rep(1e-5, parlen)))
  }
  
  likfun <- function(params){
    if(repar) params[2] <- exp(params[2])
    linERRloglikcomposite(params, data=data, doses=doses,set=set,status=status,loc=loc,in_cohort=in_cohort,corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, ccmethod=ccmethod, weights=weights)
  }
  
  
  names <- c("(Intercept)",ifelse(repar,"xi","beta"), paste0("alpha",2:(length(doses)/2)),names(data)[corrvars])
  names(initpars) <- names
  parnames(likfun) <- names
  
  
  if(ccmethod=="CCAL"){
    fxd <- NULL
    if(opt_method!="Nelder-Mead"){
      lw <- sapply(names,FUN=function(x) -Inf)
      if(!repar) lw[2] <- -1/max(data[,doses])+.0001
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else {
      lw <- NULL
      up <- NULL
    }
  } else if (ccmethod=="CL"){
    fxd <- sapply(names(data)[corrvars], function(i) 0)
    if(opt_method!="Nelder-Mead"){
      lw <- sapply(names[!names %in% names(data)[corrvars]],FUN=function(x) -Inf)
      if(!repar) lw[2] <- -1/max(data[data[,status]==1,doses])+.0001
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else {
      lw <- NULL
      up <- NULL
    }
  } else if(ccmethod=="CCML"){
    if(opt_method=="L-BFGS-B"){
      fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
      lw <- sapply(c(ifelse(repar,"xi","beta"),tail(names, length(corrvars))),FUN=function(x) -Inf)
      if(!repar) {lw[1] <- -1/max(data[,doses][cbind(1:nrow(data), data[,loc])])+.0001}
      #if(length(corrvars)==0){
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else { # Nelder-Mead
      fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
      lw <-NULL
      up <- NULL
    }
    #} else{
    #  up <- NULL
    #}
  } else if(ccmethod=="meandose"){
    if(opt_method=="L-BFGS-B"){
      fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
      lw <- sapply(c(ifelse(repar,"xi","beta"),tail(names, length(corrvars))),FUN=function(x) -Inf)
      if(!repar) {lw[1] <- -1/max(rowMeans(data[,doses]))+.0001}
      #if(length(corrvars)==0){
      if(repar){
        up <- list(xi=log(uplimBeta))
      } else {
        up <- list(beta=uplimBeta)
      }
    } else { # Nelder-Mead
      fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
      up <- NULL
      lw <- NULL
    }
    #} else{
    #  up <- NULL
    #}
  }
  
  
  fit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd,upper=up, lower=lw, control=fitopt, method=opt_method)
  
  if(fitNull){ # This needs to be checked
    if(opt_method!="Brent"){
      fitoptnull <- fitopt
      fitoptnull$ndeps <- fitoptnull$ndeps[-2]
      if(repar){
        fxd2 <- c(fxd, list(xi=-999999))
      } else {
        fxd2 <- c(fxd, list(beta=0))
      }
    } else {
      fitoptnull <- fitopt
      fxd2 <- list(params=ifelse(repar, -999999, 0))
    }
    
    lw2 <- lw[-1]
    up2 <- up[-1]
    if(length(lw2)==0) lw2 <- NULL
    if(length(up2)==0) up2 <- NULL
    nullfit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd2,upper=up2, lower=lw2, control=fitoptnull, method=opt_method)
  } else {
    nullfit <- NULL
  }
  
  proflik <- function(x){
    if(opt_method!="Brent"){
      fitoptprof <- fitopt
      fitoptprof$ndeps <- fitoptprof$ndeps[-1]
      if(repar){
        fxd3 <- c(fxd, list(xi=x))
      } else {
        fxd3 <- c(fxd, list(beta=x))
      }
    } else {
      fitoptprof <- fitopt
      fxd3 <- list(params=x)
    }
    lw3 <- lw[-1]
    up3 <- up[-1]
    if(length(lw3)==0) lw3 <- NULL
    if(length(up3)==0) up3 <- NULL
    
    mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd3,upper=up3, lower=lw3, control=fitoptprof, method=opt_method)@min
  }
  
  
  list(fit=fit, nullfit=nullfit, proflik=proflik)
  
}


linearERRfitcompositeexp <- function(data, doses, set, status, loc, in_cohort, corrvars=NULL,corrvars_exclude_conditional=NULL, ccmethod="CCAL", initpars=rep(0,1+length(doses)/2+length(corrvars)), fitopt=list(maxit=5000), fitNull=TRUE, useOld=FALSE, uplimBeta=5, weights){
  
  if(ccmethod=="CL") corrvars <- NULL
  
  opt_method <- "Nelder-Mead"
  
  
  if(is.null(fitopt$maxit)){
    fitopt <- c(fitopt, list(maxit=5000))
  }
  if(is.null(fitopt$reltol) & opt_method!="L-BFGS-B"){
    fitopt <- c(fitopt, list(reltol=1e-10))
  }
  if (is.null(fitopt$pgtol) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(pgtol=1e-8))
  }
  if (is.null(fitopt$factr) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(factr=1e4))
  }
  # if (is.null(fitopt$ndeps) & opt_method=="L-BFGS-B"){
  #   parlen <- 2+ifelse(ccmethod %in% c("CCAL","CL"),length(doses)/2-1,0)+length(corrvars)
  #   fitopt <- c(fitopt, list(ndeps=rep(1e-5, parlen)))
  # }
  
  likfun <- function(params){
    #if(repar) params[2] <- exp(params[2])
    linERRloglikcompositeexp(params, data=data, doses=doses,set=set,status=status,loc=loc,in_cohort=in_cohort,corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, ccmethod=ccmethod, weights=weights)
  }
  
  
  names <- c("(Intercept)","beta", paste0("alpha",2:(length(doses)/2)),names(data)[corrvars])
  names(initpars) <- names
  parnames(likfun) <- names
  
  
  if(ccmethod=="CCAL"){
    fxd <- NULL
    lw <- NULL
    up <- NULL
    
  } else if (ccmethod=="CL"){
    fxd <- sapply(names(data)[corrvars], function(i) 0)
    lw <- NULL
    up <- NULL
    
  } else if(ccmethod=="CCML"){
    
    fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
    lw <-NULL
    up <- NULL
    
  } else if(ccmethod=="meandose"){
    fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
    up <- NULL
    lw <- NULL
    
  }
  
  
  fit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd,upper=up, lower=lw, control=fitopt, method=opt_method)
  
  if(fitNull){ # This needs to be checked
    if(opt_method!="Brent"){
      fitoptnull <- fitopt
      fitoptnull$ndeps <- fitoptnull$ndeps[-2]
      
      fxd2 <- c(fxd, list(beta=0))
      
    } 
    
    lw2 <- lw[-1]
    up2 <- up[-1]
    if(length(lw2)==0) lw2 <- NULL
    if(length(up2)==0) up2 <- NULL
    nullfit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd2,upper=up2, lower=lw2, control=fitoptnull, method=opt_method)
  } else {
    nullfit <- NULL
  }
  
  proflik <- function(x){
    if(opt_method!="Brent"){
      fitoptprof <- fitopt
      fitoptprof$ndeps <- fitoptprof$ndeps[-1]
      
      fxd3 <- c(fxd, list(beta=x))
      
    } else {
      fitoptprof <- fitopt
      fxd3 <- list(params=x)
    }
    lw3 <- lw[-1]
    up3 <- up[-1]
    if(length(lw3)==0) lw3 <- NULL
    if(length(up3)==0) up3 <- NULL
    
    mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd3,upper=up3, lower=lw3, control=fitoptprof, method=opt_method)@min
  }
  
  
  list(fit=fit, nullfit=nullfit, proflik=proflik)
  
}


linearERRfitcompositecat <- function(data, doses, set, status, loc, in_cohort, corrvars=NULL,corrvars_exclude_conditional=NULL, ccmethod="CCAL", initpars=rep(0,1+length(doses)/2+length(corrvars)), fitopt=list(maxit=5000), fitNull=TRUE, useOld=FALSE, uplimBeta=5, weights){
  
  if(ccmethod=="CL") corrvars <- NULL
  
  opt_method <- "Nelder-Mead"
  
  
  if(is.null(fitopt$maxit)){
    fitopt <- c(fitopt, list(maxit=5000))
  }
  if(is.null(fitopt$reltol) & opt_method!="L-BFGS-B"){
    fitopt <- c(fitopt, list(reltol=1e-10))
  }
  if (is.null(fitopt$pgtol) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(pgtol=1e-8))
  }
  if (is.null(fitopt$factr) & opt_method=="L-BFGS-B"){
    fitopt <- c(fitopt, list(factr=1e4))
  }
  # if (is.null(fitopt$ndeps) & opt_method=="L-BFGS-B"){
  #   parlen <- 2+ifelse(ccmethod %in% c("CCAL","CL"),length(doses)/2-1,0)+length(corrvars)
  #   fitopt <- c(fitopt, list(ndeps=rep(1e-5, parlen)))
  # }
  
  likfun <- function(params){
    #if(repar) params[2] <- exp(params[2])
    linERRloglikcompositecat(params, data=data, doses=doses,set=set,status=status,loc=loc,in_cohort=in_cohort,corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, ccmethod=ccmethod, weights=weights)
  }
  
  
  names <- c("(Intercept)","beta1","beta2", paste0("alpha",2:(length(doses)/2)),names(data)[corrvars])
  names(initpars) <- names
  parnames(likfun) <- names
  
  
  if(ccmethod=="CCAL"){
    fxd <- NULL
    lw <- NULL
    up <- NULL
    
  } else if (ccmethod=="CL"){
    fxd <- sapply(names(data)[corrvars], function(i) 0)
    lw <- NULL
    up <- NULL
    
  } else if(ccmethod=="CCML"){
    
    fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
    lw <-NULL
    up <- NULL
    
  } else if(ccmethod=="meandose"){
    fxd <- sapply(paste0("alpha",2:(length(doses)/2)), function(i) 0)
    up <- NULL
    lw <- NULL
    
  }
  
  
  fit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd,upper=up, lower=lw, control=fitopt, method=opt_method)
  
  if(fitNull){ # This needs to be checked
    if(opt_method!="Brent"){
      fitoptnull <- fitopt
      fitoptnull$ndeps <- fitoptnull$ndeps[-2]
      
      fxd2 <- c(fxd, list(beta1=0, beta2=0))
      
    } 
    
    lw2 <- lw[-1]
    up2 <- up[-1]
    if(length(lw2)==0) lw2 <- NULL
    if(length(up2)==0) up2 <- NULL
    nullfit <- mle2(likfun, start=initpars, vecpar=(opt_method!="Brent"),fixed=fxd2,upper=up2, lower=lw2, control=fitoptnull, method=opt_method)
  } else {
    nullfit <- NULL
  }
  
  proflik <- NULL
  
  
  list(fit=fit, nullfit=nullfit, proflik=proflik)
  
}


linearERRcomposite <- function(data, doses, set, status, loc, in_cohort, corrvars=NULL, corrvars_exclude_conditional=NULL, ccmethod="CCAL",repar=FALSE, initpars=rep(0,1+length(doses)/2+length(corrvars)), fitopt=NULL, uplimBeta=5,profCI=TRUE,doJK1=FALSE,doJK2=FALSE,jkscorethresh=.01,jkvalrange=c(-Inf, Inf), weights){
  
  if(doJK2) doJK1 <- TRUE
  
  mainfit <- linearERRfitcomposite(data=data, doses=doses, set=set,status=status, loc=loc,in_cohort=in_cohort, corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, repar=repar, initpars=initpars, fitopt=fitopt, ccmethod=ccmethod, uplimBeta=uplimBeta, weights=weights, fitNull = FALSE)
  
  likfun <- function(params){
    if(repar) params[2] <- exp(params[2])
    linERRloglikcomposite(params, data=data, doses=doses,set=set,status=status,loc=loc,in_cohort=in_cohort,corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, ccmethod=ccmethod, weights=weights)
  }
  
  MLEscore <- grad(likfun, mainfit$fit@coef, method.args=list(r=10))
  
  
  
  
  
  pval <- 1#pchisq(2*(mainfit$nullfit@min-mainfit$fit@min), df=1, lower.tail=FALSE)
  
  
  if(profCI){
    
    
    g <- function(para){
      1-pchisq(2*(mainfit$proflik(para)-mainfit$fit@min),df=1)-.05
    }
    lowLim <- tryCatch(uniroot(g, lower=ifelse(repar,-20,mainfit$fit@call$lower[1]), upper=mainfit$fit@coef[1], extendInt="no")$root, error=function(e) NA)
    upLim <- tryCatch(uniroot(g, lower=mainfit$fit@coef[1],upper=ifelse(repar,log(100),100), extendInt="no", maxiter=150)$root, error=function(e) NA)
    
  } else{
    lowLim <- NULL
    upLim <- NULL
  }
  
  MLE <- list(coef=mainfit$fit@coef,sd=sqrt(diag(mainfit$fit@vcov)), vcov=mainfit$fit@vcov, score=MLEscore, convergence=mainfit$fit@details$convergence, message=mainfit$fit@details$message, dosepval=pval, profCI=c(lo=lowLim, up=upLim), fitobj=mainfit)
  
  
  jackknife <- NULL
  
  list(MLE=MLE, jackknife=jackknife)
  
}


linearERRcompositeexp <- function(data, doses, set, status, loc, in_cohort, corrvars=NULL, corrvars_exclude_conditional=NULL, ccmethod="CCAL", initpars=rep(0,1+length(doses)/2+length(corrvars)), fitopt=NULL, uplimBeta=5,profCI=TRUE,doJK1=FALSE,doJK2=FALSE,jkscorethresh=.01,jkvalrange=c(-Inf, Inf), weights){
  
  if(doJK2) doJK1 <- TRUE
  
  mainfit <- linearERRfitcompositeexp(data=data, doses=doses, set=set,status=status, loc=loc,in_cohort=in_cohort, corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, initpars=initpars, fitopt=fitopt, ccmethod=ccmethod, uplimBeta=uplimBeta, weights=weights, fitNull = FALSE)
  
  likfun <- function(params){
    #if(repar) params[2] <- exp(params[2])
    linERRloglikcompositeexp(params, data=data, doses=doses,set=set,status=status,loc=loc,in_cohort=in_cohort,corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, ccmethod=ccmethod, weights=weights)
  }
  
  MLEscore <- grad(likfun, mainfit$fit@coef, method.args=list(r=10))
  
  
  
 
  
  pval <- 1#pchisq(2*(mainfit$nullfit@min-mainfit$fit@min), df=1, lower.tail=FALSE)
  
  
  if(profCI){
    
    
    g <- function(para){
      1-pchisq(2*(mainfit$proflik(para)-mainfit$fit@min),df=1)-.05
    }
    lowLim <- tryCatch(uniroot(g, lower=ifelse(repar,-20,mainfit$fit@call$lower[1]), upper=mainfit$fit@coef[1], extendInt="no")$root, error=function(e) NA)
    upLim <- tryCatch(uniroot(g, lower=mainfit$fit@coef[1],upper=ifelse(repar,log(100),100), extendInt="no", maxiter=150)$root, error=function(e) NA)
    
  } else{
    lowLim <- NULL
    upLim <- NULL
  }
  
  MLE <- list(coef=mainfit$fit@coef,sd=sqrt(diag(mainfit$fit@vcov)), vcov=mainfit$fit@vcov, score=MLEscore, convergence=mainfit$fit@details$convergence, message=mainfit$fit@details$message, dosepval=pval, profCI=c(lo=lowLim, up=upLim), fitobj=mainfit)
  
  
  jackknife <- NULL
  
  list(MLE=MLE, jackknife=jackknife)
  
}


linearERRcompositecat <- function(data, doses, set, status, loc, in_cohort, corrvars=NULL, corrvars_exclude_conditional=NULL, ccmethod="CCAL", initpars=rep(0,2+length(doses)/2+length(corrvars)), fitopt=NULL, uplimBeta=5,profCI=TRUE,doJK1=FALSE,doJK2=FALSE,jkscorethresh=.01,jkvalrange=c(-Inf, Inf), weights){
  
  if(doJK2) doJK1 <- TRUE
  
  mainfit <- linearERRfitcompositecat(data=data, doses=doses, set=set,status=status, loc=loc,in_cohort=in_cohort, corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, initpars=initpars, fitopt=fitopt, ccmethod=ccmethod, uplimBeta=uplimBeta, weights=weights, fitNull = FALSE)
  
  likfun <- function(params){
    #if(repar) params[2] <- exp(params[2])
    linERRloglikcompositecat(params, data=data, doses=doses,set=set,status=status,loc=loc,in_cohort=in_cohort,corrvars=corrvars,corrvars_exclude_conditional=corrvars_exclude_conditional, ccmethod=ccmethod, weights=weights)
  }
  
  MLEscore <- grad(likfun, mainfit$fit@coef, method.args=list(r=10))
  
  
  
  
  
  pval <- 1#pchisq(2*(mainfit$nullfit@min-mainfit$fit@min), df=1, lower.tail=FALSE)
  
  
  if(profCI){
    
    
    g <- function(para){
      1-pchisq(2*(mainfit$proflik(para)-mainfit$fit@min),df=1)-.05
    }
    lowLim <- tryCatch(uniroot(g, lower=ifelse(repar,-20,mainfit$fit@call$lower[1]), upper=mainfit$fit@coef[1], extendInt="no")$root, error=function(e) NA)
    upLim <- tryCatch(uniroot(g, lower=mainfit$fit@coef[1],upper=ifelse(repar,log(100),100), extendInt="no", maxiter=150)$root, error=function(e) NA)
    
  } else{
    lowLim <- NULL
    upLim <- NULL
  }
  
  MLE <- list(coef=mainfit$fit@coef,sd=sqrt(diag(mainfit$fit@vcov)), vcov=mainfit$fit@vcov, score=MLEscore, convergence=mainfit$fit@details$convergence, message=mainfit$fit@details$message, dosepval=pval, profCI=c(lo=lowLim, up=upLim), fitobj=mainfit)
  
  
  jackknife <- NULL
  
  list(MLE=MLE, jackknife=jackknife)
  
}

















