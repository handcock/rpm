#' @keywords internal
#' @importFrom MASS ginv
rpm.hessian <- function(theta,Sd,Xd,Zd,NumBeta,NumGamma,NumGammaW,NumGammaM,pmfW,pmfM,pmf,counts,gw,gm,N,sampling,constraints,verbose=TRUE){
   if(length(N) > 1){
     ext.covar <- diag(length(theta))
     dimnames(ext.covar) <- list(names(theta),names(theta))
     return(list(covar=ext.covar[1:NumBeta,1:NumBeta],ext.covar=ext.covar,covar.unconstrained=ext.covar))
   }
   hloglikfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, sampling, constraints){
      NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
      beta <- theta[1:NumBeta]
      GammaW <- theta[NumBeta+(1:NumGammaW)]
      GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
      gw <- theta[(NumBeta+NumGammaW+NumGammaM+1)]
      gm <- log(1-exp(gw))
      -hloglik(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
   H <- hloglikfun(theta,
     Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
     pmfW=pmfW, pmfM=pmfM, pmf=pmf, counts=counts, gw=gw, gm=gm, sampling=sampling, constraints=constraints)
#
   jeqfun <- function(theta, Sd, Xd, Zd, NumGammaW, NumGammaM, pmfW, pmfM, pmf, counts, gw, gm, sampling, constraints){
     NumBeta <- dim(Sd)[3]+dim(Xd)[3]+dim(Zd)[3]
     beta <- theta[1:NumBeta]
     GammaW <- theta[NumBeta+(1:NumGammaW)]
     GammaM <- theta[(NumBeta+NumGammaW)+(1:NumGammaM)]
     gw <- theta[(NumBeta+NumGammaW+NumGammaM+1)]
     gm <- log(1-exp(gw))
     geqcond(beta, GammaW, GammaM, Sd, Xd, Zd, dim(Sd), dim(Xd), dim(Zd), pmfW, pmfM, pmf, counts, gw, gm, constraints)
    }
   GJ <- jeqfun_default(theta,
           Sd=Sd,Xd=Xd,Zd=Zd,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
           pmfW=pmfW, pmfM=pmfM, pmf=pmf, counts=counts, gw=gw, gm=gm, sampling=sampling,
           constraints=constraints)
   dimnames(H) <- list(names(theta),names(theta))
   if(constraints==1){
     dimnames(GJ) <- list(c(names(theta)[(NumBeta+1):(NumBeta+NumGamma)]),
                          names(theta))
   }else{
     dimnames(GJ) <- list(c(names(theta)[(NumBeta+1):(NumBeta+NumGamma)],paste0("M_",1:(NumGammaM))),
                          names(theta))
   }
   dimnames(H) <- list(names(theta)[1:(NumBeta+NumGamma+1)], names(theta)[1:(NumBeta+NumGamma+1)])
   Hi <- try(MASS::ginv(H))
   dimnames(Hi) <- list(names(theta)[1:(NumBeta+NumGamma+1)], names(theta)[1:(NumBeta+NumGamma+1)])
   if(inherits(Hi,"try-error")){
     if(verbose) message("Trouble computing the standard errors. They are approximate.")
     Hi <- diag(1/diag(H))
   }
   # See Hartmann and Hartwig (1996)
   Md <- rbind(H,GJ)
   Md <- cbind(Md,rbind(t(GJ),matrix(0,ncol=nrow(GJ),nrow=nrow(GJ))))
   Mi=try(MASS::ginv(Md))
   if(inherits(Mi,"try-error")){
     if(verbose) message("Trouble computing the standard errors. They are approximate.")
     Mi <- diag(1/diag(Md))
   }
   V=Mi[1:NumBeta,1:NumBeta]
   a <- diag(V)
   a[is.na(a)] <- 0
   if(all(is.na(diag(V)) | abs(a)<1e-15 | abs(a) > 100)){
     if(verbose) message("Trouble computing the standard errors. Using the unconstrained versions.")
     covar <- Hi[1:NumBeta,1:NumBeta]
   }else{
     covar <- V
   }
# Next constrained
  ext.covar=Mi[1:(NumBeta+NumGamma+1),1:(NumBeta+NumGamma+1)]
  dimnames(ext.covar) <- list(names(theta),names(theta))
  list(covar=covar,ext.covar=ext.covar,covar.unconstrained=Hi)
}

#' @keywords internal
.catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 
#' @keywords internal
mode.density <- function(x){
   minx <- min(x)
   maxx <- max(x)
   xp <- seq(minx, maxx, length=1000)
   a=bgk_kde(x,n=2^(ceiling(log((maxx-minx))/log(2))),MIN=minx,MAX=maxx)
   posdens <- stats::spline(x=a[1,],y=a[2,],xout=xp)$y
   posdens <- 1000*posdens / ((maxx-minx)*sum(posdens))
   map <- xp[which.max(posdens)]
   map
}

#' @keywords internal
is.empty <- function(x, mode=NULL){
    if (is.null(x)) return(TRUE)
    if (is.null(mode)) mode <- class(x)
    identical(vector(mode,1),c(x,vector(class(x),1)))
}
