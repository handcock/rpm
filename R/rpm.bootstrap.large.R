rpm.bootstrap.large <- function(I, solution,
                    S,X,Z,sampling_design,Xdata,Zdata,sampled,Xid,Zid,pair_id,X_w,Z_w,
                    num_sampled,
                    NumBeta,NumGammaW,NumGammaM,num_Xu,num_Zu,cnW, cnM,LB,UB,
                    control){

      X_w_rel <- rep(1,nrow(Xdata))
      Z_w_rel <- rep(1,nrow(Zdata))
      if(sampling_design == "stock-stock"){
        paired_W <- !is.na(Xdata[,pair_id]) 
        M_paired_to_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
        I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
        X_w_rel[paired_W] <- 0.0
        X_w_rel[paired_W][ I] <- 1
        Z_w_rel[M_paired_to_W] <- 0.0
        Z_w_rel[M_paired_to_W][!I] <- 1
        I <- sample(c(Xdata[,Xid],
                      Zdata[,Zid]),
                      prob=c(X_w_rel,Z_w_rel),
                    replace=TRUE,size=num_sampled)
      }
      if(sampling_design == "census"){
        paired_W <- !is.na(Xdata[,pair_id]) 
        M_paired_to_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
        I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
        I <- c(Xdata[paired_W,Xid][I],Zdata[M_paired_to_W,Zid][!I])
        I <- c(I,sample(c(Xdata[!paired_W,Xid], Zdata[is.na(Zdata[,pair_id]),Zid]),
             replace=TRUE,size=num_sampled-sum(paired_W)))
      }
      if(sampling_design == "stock-flow"){
        paired_W <- !is.na(Xdata[,pair_id]) 
        M_paired_to_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
        I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
        X_w_rel[paired_W] <- 0.0
        X_w_rel[paired_W][ I] <- 2
        Z_w_rel[M_paired_to_W] <- 0.0
        Z_w_rel[M_paired_to_W][!I] <- 2
        I <- sample(c(Xdata[,Xid],
                      Zdata[,Zid]),
                      prob=c(X_w_rel,Z_w_rel),
                    replace=TRUE,size=num_sampled)
      }
      Xmatch <- match(I,Xdata[,Xid])
      Xmatch <- Xmatch[!is.na(Xmatch)]
      Zmatch <- match(I,Zdata[,Zid])
      Zmatch <- Zmatch[!is.na(Zmatch)]
      XdataS <- Xdata[Xmatch[!is.na(Xmatch)],]
      ZdataS <- Zdata[Zmatch[!is.na(Zmatch)],]
  #   Find the people paired to the sampled people
      paired_and_sampled_W <- !is.na(XdataS[,pair_id])
      M_paired_to_sampled_W <- match(XdataS[paired_and_sampled_W,pair_id], Zdata[,Zid])
      paired_and_sampled_M <- !is.na(ZdataS[,pair_id])
      W_paired_to_sampled_M <- match(ZdataS[paired_and_sampled_M,pair_id], Xdata[,Xid])
      XdataM <- Xdata[W_paired_to_sampled_M,]
      ZdataW <- Zdata[M_paired_to_sampled_W,]
  
      if(sampling_design == "stock-flow"){
        if(nrow(XdataS)>0) XdataS[,sampled] <- TRUE
        if(nrow(ZdataS)>0) ZdataS[,sampled] <- TRUE
        if(nrow(XdataS)>0) XdataS[paired_and_sampled_W,X_w] <- XdataS[paired_and_sampled_W,X_w] + ZdataW[,Z_w]
        if(nrow(ZdataS)>0) ZdataS[paired_and_sampled_M,Z_w] <- ZdataS[paired_and_sampled_M,Z_w] + XdataM[,X_w]
        if(nrow(XdataM)>0) XdataM[,sampled] <- FALSE
        if(nrow(ZdataW)>0) ZdataW[,sampled] <- FALSE
      } else if (sampling_design == "stock-stock"){ 
        if(nrow(XdataS)>0) XdataS[,sampled] <- TRUE
        if(nrow(ZdataS)>0) ZdataS[,sampled] <- TRUE
        if(nrow(XdataM)>0) XdataM[,sampled] <- TRUE
        if(nrow(ZdataW)>0) ZdataW[,sampled] <- TRUE
      }

      Xdata=rbind(XdataS,XdataM)
      Zdata=rbind(ZdataS,ZdataW)

      # IDs of the women matched to the sampled men (and vice versa)
      if(sampling_design != "census"){
        paired_and_sampled_M <- Zdata[,sampled] & !is.na(Zdata[,pair_id])
        paired_and_sampled_W <- Xdata[,sampled] & !is.na(Xdata[,pair_id])
      } else {
	paired_and_sampled_M <- !is.na(Zdata[,pair_id])
        paired_and_sampled_W <- !is.na(Xdata[,pair_id])    
      }
      M_paired_to_sampled_W <- match(Xdata[paired_and_sampled_W,pair_id], Zdata[,Zid])
      W_paired_to_sampled_M <- match(Zdata[paired_and_sampled_M,pair_id], Xdata[,Xid])
      XdataM <- Xdata[W_paired_to_sampled_M,]
      ZdataW <- Zdata[M_paired_to_sampled_W,]

      if (sampling_design == "census") {
       pmfW = as.numeric(stats::xtabs(~ factor(Xtype,1:num_Xu), data=Xdata))
       pmfM = as.numeric(stats::xtabs(~ factor(Ztype,1:num_Zu), data=Zdata))
      }
      if (sampling_design == "stock-stock") {
       subset=Xdata[,sampled] &  is.na(Xdata[,pair_id])
       pmfW_S = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=Xdata, subset=subset))
       subset=Xdata[,sampled] & !is.na(Xdata[,pair_id])
       pmfW_P = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=Xdata, subset=subset))
       pmfW = pmfW_S + pmfW_P
       subset=Zdata[,sampled] &  is.na(Zdata[,pair_id])
       pmfM_S = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=Zdata, subset=subset))
       subset=Zdata[,sampled] & !is.na(Zdata[,pair_id])
       pmfM_P = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=Zdata, subset=subset))
       pmfM = pmfM_S + pmfM_P
      }
      if (sampling_design == "stock-flow") {
        pmfW = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=Xdata, subset=sampled))
        subset=Zdata[,sampled] & !is.na(Zdata[,pair_id])
        pmfW = pmfW + as.numeric(stats::xtabs(Zdata$Z_w[subset] ~ factor(Xdata$Xtype[W_paired_to_sampled_M],1:num_Xu)))
        pmfM = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=Zdata, subset=sampled))
        subset=Xdata[,sampled] & !is.na(Xdata[,pair_id])
        pmfM = pmfM + as.numeric(stats::xtabs(Xdata$X_w[subset] ~ factor(Zdata$Ztype[M_paired_to_sampled_W],1:num_Zu)))
      }
      pmfW = pmfW/sum(pmfW)
      pmfM = pmfM/sum(pmfM)
      names(pmfW) <- cnW
      names(pmfM) <- cnM
  
      if (sampling_design == "census") {
        X_sel <- is.na(Xdata[,pair_id])
        Z_sel <- is.na(Zdata[,pair_id])
      }else{
        X_sel <- Xdata[,sampled] & is.na(Xdata[,pair_id])
        Z_sel <- Zdata[,sampled] & is.na(Zdata[,pair_id])
      }
      Xcounts_single = as.numeric(stats::xtabs(~ factor(Xtype, 1:num_Xu), data=Xdata, subset=X_sel)) 
      Zcounts_single = as.numeric(stats::xtabs(~ factor(Ztype, 1:num_Zu), data=Zdata, subset=Z_sel))
  
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
      
      # The number of people in the population
      if(sampling_design == "stock-stock"){
        N_w = sum(Xdata[Xdata[,sampled] & is.na(Xdata[,pair_id]),X_w])
        N_w = N_w + sum(Zdata[Zdata[,sampled] & !is.na(Zdata[,pair_id]),Z_w]) # The population size
        N_m = sum(Zdata[Zdata[,sampled] & is.na(Zdata[,pair_id]),Z_w])
        N_m = N_m + sum(Xdata[Xdata[,sampled] & !is.na(Xdata[,pair_id]),X_w]) # The population size
        N = N_w + N_m
        gw = log(N_w/N) # to ensure exp(gw)+exp(gm) = 1
        gm = log(N_m/N) # to ensure exp(gw)+exp(gm) = 1
      }
      if(sampling_design == "stock-flow"){
        N_w = sum(Xdata[Xdata[,sampled],X_w])
        N_m = sum(Zdata[Zdata[,sampled],Z_w])
        N = N_w + N_m
        gw = log(N_w/N) # to ensure exp(gw)+exp(gm) = 1
        gm = log(N_m/N) # to ensure exp(gw)+exp(gm) = 1
      }
      if(sampling_design == "census"){
        N = nrow(Xdata) + nrow(Zdata) # The population size
        gw = log(nrow(Xdata)/N) # to ensure exp(gw)+exp(gm) = 1
        gm = log(nrow(Zdata)/N)
      }

      # Compute the number sampled (mainly for the s.e. computation)
      if (sampling_design == "stock-stock") {
       num_sampled <- nrow(Xdata)+nrow(Zdata)-sum(!is.na(Zdata[,pair_id]))
      }else{
       if (sampling_design == "stock-flow") {
        num_sampled <- sum(Xdata[,sampled])+sum(Zdata[,sampled])
       }else{
        # census for now
        num_sampled <- nrow(Xdata)+nrow(Zdata)
       }
      }

      pmf = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(pmf) <- c(cnM,"singles")
      rownames(pmf) <- c(cnW,"singles")
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
  
      if (sampling_design == "census") {
        Xtype_single = as.numeric(stats::xtabs(~ factor(Xtype, 1:num_Xu), data=Xdata, subset=X_sel)) 
        Ztype_single = as.numeric(stats::xtabs(~ factor(Ztype, 1:num_Zu), data=Zdata, subset=Z_sel))
        pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu) + factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu) + factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
        pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N
        counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      }
      if (sampling_design == "stock-stock") {
        Xtype_single = as.numeric(stats::xtabs(X_w ~ factor(Xtype, 1:num_Xu), data=Xdata, subset=X_sel)) 
        Ztype_single = as.numeric(stats::xtabs(Z_w ~ factor(Ztype, 1:num_Zu), data=Zdata, subset=Z_sel))
        pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(Zdata[M_paired_to_sampled_W,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(Xdata[W_paired_to_sampled_M,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
        pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N
        counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(Zdata[M_paired_to_sampled_W,"Ztype"],1:num_Zu)))
        counts[1:num_Xu,1:num_Zu] <- 2*counts[1:num_Xu,1:num_Zu]
      }
      if (sampling_design == "stock-flow") {
        Xtype_single = as.numeric(stats::xtabs(X_w ~ factor(Xtype, 1:num_Xu), data=Xdata, subset=X_sel)) 
        Ztype_single = as.numeric(stats::xtabs(Z_w ~ factor(Ztype, 1:num_Zu), data=Zdata, subset=Z_sel))
        pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(Zdata[M_paired_to_sampled_W,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(Xdata[W_paired_to_sampled_M,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
        pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N
        counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(Zdata[M_paired_to_sampled_W,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(Xdata[W_paired_to_sampled_M,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      }

      if (!is.empty(Xcounts_single)) {
        pmf[1:num_Xu,1+num_Zu] = Xtype_single / N
        counts[1:num_Xu,1+num_Zu] = Xcounts_single
      }
      if (!is.empty(Zcounts_single)) {
        pmf[1+num_Xu,1:num_Zu] = Ztype_single / N
        counts[1+num_Xu,1:num_Zu] = Zcounts_single
      }
  
      if (sampling_design == "census") { 
        pmfN <- pmf*num_sampled
      }else{
        pmfN <- counts
      }

      control$xtol_rel=control$bs.xtol_rel
      control$maxeval=control$bs.maxeval
      
      out.text <- capture.output(
       out.fit <- nloptr::nloptr(x0=solution, eval_f=loglikfun_nog, 
                 eval_grad_f=gloglikfun_nog,
                 eval_g_eq=eqfun_nog,  eval_jac_g_eq=jeqfun_nog,
                 lb=LB,ub=UB,
                 Sd=S,Xd=X,Zd=Z,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                 pmfW=pmfW, pmfM=pmfM, pmf=pmf, counts=pmfN, gw=gw, gm=gm, N=N,
                 sampling=sampling_design, constraints=1,
                 opts=control)
      )
      if(any(startsWith("Error",c(" ",out.text)))){
        message(sprintf("Optimization for starting at the MLPLE is overly constrained. Estimates may be unstable."))
      }

      # Remove the large associated environments
      out.fit$nloptr_environment<- NULL
      out.fit$eval_g_eq <- NULL
      out.fit$eval_f  <- NULL

      th_hat <-  out.fit$solution
      pmf_est <- exp(augpmfnew(th_hat[1:NumBeta],
                GammaW=th_hat[NumBeta+(1:NumGammaW)], 
                GammaM=th_hat[(NumBeta+NumGammaW)+(1:NumGammaM)],
                S, X, Z,
                pmfW, pmfM, gw=gw, gm=gm))
      pmf_est[nrow(pmf_est),ncol(pmf_est)] <- 0
      pmf_est[-nrow(pmf_est), -ncol(pmf_est)] <- 2*pmf_est[-nrow(pmf_est), -ncol(pmf_est)]
      pmf_est <- pmf_est/sum(pmf_est)

      PMF_SW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
      PMF_SW <- PMF_SW / (PMF_SW + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],1,sum))
      PMF_SM <- pmf_est[ nrow(pmf_est),-ncol(pmf_est),drop=FALSE]
      PMF_SM <- PMF_SM / (PMF_SM + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],2,sum))
      PMF_PW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
      PMF_PW <- 1 - PMF_SW
      PMF_PM <- 1 - PMF_SM 

      names(PMF_SW) <- paste0("PMF_Single.W.",cnW)
      names(PMF_SM) <- paste0("PMF_Single.M.",cnM)
      names(PMF_PW) <- paste0("PMF_Partnered.W.",cnW)
      names(PMF_PM) <- paste0("PMF_Partnered.M.",cnM)
      if(control$logodds_single){
        b <- sum(PMF_SW*pmfW)
        LOGODDS_SW <- th_hat[(NumBeta+1):(NumBeta+NumGammaW)] - log(b/(1-b))
        b <- sum(PMF_SM*pmfM)
        LOGODDS_SM <- th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGammaW+NumGammaM)] - log(b/(1-b))
      }else{
        LOGODDS_SW <- log(PMF_SW/(1-PMF_SW)) #+ log(2)  # VIP note the log(2) hack. 
        LOGODDS_SM <- log(PMF_SM/(1-PMF_SM)) #+ log(2)
      }
      names(LOGODDS_SW) <- paste0("LOD_Single.W.",cnW)
      names(LOGODDS_SM) <- paste0("LOD_Single.M.",cnM)

      pmfN_households <- pmfN
      pmfN_households[-nrow(pmfN), -ncol(pmfN)] <- 0.5*pmfN[-nrow(pmfN), -ncol(pmfN)]
      pmf_households <- pmfN_households / sum(pmfN_households)
      loglik <- stats::dmultinom(x=pmfN_households,prob=pmf_est,log=TRUE)

      list(est=th_hat, LOGODDS_SW=LOGODDS_SW,LOGODDS_SM=LOGODDS_SM,
           Xdata=Xdata,
           Zdata=Zdata,
           pmfW=pmfW, pmfM=pmfM,
           pmf=pmf, counts=counts, nobs=num_sampled,
           pmf_est=pmf_est,
           aic = 2*NumBeta-2*loglik, bic=log(num_sampled)*NumBeta-2*loglik, loglik=loglik
          )
     }
