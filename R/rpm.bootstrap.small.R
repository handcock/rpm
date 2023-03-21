rpm.bootstrap.small <- function(i, solution, num_women, num_men, Jw, Jm, U_star, V_star,
     S, X, Z, pmfW, pmfM,
     Xu, Zu, num_Xu, num_Zu, cnW, cnM, Xid, Zid, pair_id, sampled, sampling_design,
     NumBeta,NumGammaW,NumGammaM,LB,UB,control,
     num_sampled){

        Ws <- stats::rmultinom(n=1, size=num_women, prob=pmfW)
        Ms <- stats::rmultinom(n=1, size=num_men, prob=pmfM)
        Ws <- rep(seq_along(pmfW),Ws)
        Ms <- rep(seq_along(pmfM),Ms)
        # create utility matrices
        U_star = matrix(0, nrow=length(Ws), ncol = length(Ms))
        V_star = matrix(0, nrow=length(Ms), ncol = length(Ws))
        beta_S <- solution[1:dim(S)[3]]
        for (ii in 1:dim(S)[3]) {
          U_star = U_star + S[Ws,Ms,ii] * beta_S[ii] * 0.5
        }
        if(dim(X)[3]>0){
         beta_w <- solution[dim(S)[3]+(1:dim(X)[3])]
         for (ii in 1:dim(X)[3]) {
          U_star = U_star + X[Ws,Ms,ii] * beta_w[ii]
         }
        }
        for (ii in 1:dim(S)[3]) {
          V_star = V_star + t(S[Ws,Ms,ii]) * beta_S[ii] * 0.5
        }
        if(dim(Z)[3]>0){
         beta_m <- solution[dim(S)[3]+dim(X)[3] + (1:dim(Z)[3])]
         for (ii in 1:dim(Z)[3]) {
          V_star = V_star + Z[Ms,Ws,ii] * beta_m[ii]
         }
        }

          eta  <- -log(-log(matrix(stats::runif(num_women * num_men), nrow=num_women)))
          zeta <- -log(-log(matrix(stats::runif(num_women * num_men), nrow=num_men)))
  
          eta0  <- -log(-log(stats::runif(num_women)))+log(Jw)
          zeta0 <- -log(-log(stats::runif(num_men)))+log(Jm)
    
          ############## temp to match one-to-one#################
          U <- cbind(eta0,  U_star + eta)
          V <- cbind(zeta0, V_star + zeta)
  
          ########################################################
  
          # generate the matching (W-optimal)
          # uses Menzel's GS which allows remaining single
          ############## temp to match one-to-one#################
          mu = Gale_Shapley(U, V, return.data.frame=TRUE, nmax=10000*nrow(U))
          mu=data.frame(mu, Xtype=c(Ws, Ms), sampled=rep(TRUE,nrow(mu)))
          colnames(mu)[c(1,3)] <- c(Xid, pair_id)
          Xdata <- subset(mu, gender=="F")
          Xdata <- data.frame(Xdata,as.data.frame(Xu[Ws,-1,drop=FALSE]))
          Zdata <- subset(mu, gender=="M")
          Zdata <- data.frame(Zdata,as.data.frame(Zu[Ms,-1,drop=FALSE]))
          colnames(Zdata)[match(Xid,colnames(Zdata))] <- Zid
          colnames(Zdata)[match("Xtype",colnames(Zdata))] <- "Ztype"

          X_w_rel <- rep(1,nrow(Xdata))
          Z_w_rel <- rep(1,nrow(Zdata))
          if(sampling_design %in% c("stock-stock", "stock-flow")){
          if(sampling_design %in% c("stock-stock", "census")){
            paired_W <- !is.na(Xdata[,pair_id])
            paired_M <- !is.na(Zdata[,pair_id])
            M_paired_to_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
            I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_M), replace=TRUE)
            X_w_rel[paired_W] <- 0.0
            X_w_rel[paired_W][ I] <- 1# num_sampled / (nrow(Xdata)+nrow(Zdata))
            Z_w_rel[paired_M] <- 0.0
            Z_w_rel[M_paired_to_W][!I] <- 1.0
          }
          if(sampling_design %in% c("stock-flow")){
            paired_W <- !is.na(Xdata[,pair_id])
            M_paired_to_W <- match(Xdata[paired_W,pair_id], Zdata[,Zid])
            I <- sample(rep(c(TRUE,FALSE)), size=sum(paired_W), replace=TRUE)
            X_w_rel[paired_W] <- 0.0
            X_w_rel[paired_W][ I] <- 2
            Z_w_rel[M_paired_to_W] <- 0.0
            Z_w_rel[M_paired_to_W][!I] <- 2
          }

          if(sum(c(X_w_rel,Z_w_rel) > 0) > num_sampled){
            I <- sample(c(Xdata[,Xid], Zdata[,Zid]),
                          prob=c(X_w_rel,Z_w_rel),
                        replace=FALSE,size=num_sampled)
          }else{
            I <- sample(c(Xdata[,Xid], Zdata[,Zid]),
                          prob=c(X_w_rel,Z_w_rel),
                        replace=TRUE,size=num_sampled)
          }
          Xmatch <- match(I,Xdata[,Xid])
          Xmatch <- Xmatch[!is.na(Xmatch)]
          Zmatch <- match(I,Zdata[,Zid])
          Zmatch <- Zmatch[!is.na(Zmatch)]
          XdataS <- Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataS <- Zdata[Zmatch[!is.na(Zmatch)],]
  #       Find the people paired to the sampled people
          paired_and_sampled_W <- !is.na(XdataS[,pair_id])
          M_paired_to_sampled_W <- match(XdataS[paired_and_sampled_W,pair_id], Zdata[,Zid])
          paired_and_sampled_M <- !is.na(ZdataS[,pair_id])
          W_paired_to_sampled_M <- match(ZdataS[paired_and_sampled_M,pair_id], Xdata[,Xid])
          XdataP <- Xdata[W_paired_to_sampled_M,]
          ZdataP <- Zdata[M_paired_to_sampled_W,]

          if (sampling_design == "stock-flow") {
            if(nrow(XdataP)>0) XdataP$sampled <- rep(FALSE,nrow(XdataP))
            if(nrow(ZdataP)>0) ZdataP$sampled <- rep(FALSE,nrow(ZdataP))
           }else{
            # If not "stock-flow" then all are sampled
            XdataP$sampled <- rep(TRUE,nrow(XdataP))
            ZdataP$sampled <- rep(TRUE,nrow(ZdataP))
           }

          if(nrow(XdataS)>0) XdataS$sampled <- rep(TRUE,nrow(XdataS))
          if(nrow(ZdataS)>0) ZdataS$sampled <- rep(TRUE,nrow(ZdataS))

          Xdata <- rbind(XdataP, XdataS)
          Zdata <- rbind(ZdataP, ZdataS)

          }else{
            # a census
            Xdata$sampled <- rep(TRUE,nrow(Xdata))
            Zdata$sampled <- rep(TRUE,nrow(Zdata))
          }
  
          X_w <- "X_w"
          Z_w <- "Z_w"
          Xdata$X_w <- rep(0, nrow(Xdata))
          Xdata$X_w[Xdata$sampled] <- num_women/sum(Xdata$sampled)
          Zdata$Z_w <- rep(0, nrow(Zdata))
          Zdata$Z_w[Zdata$sampled] <- num_men/sum(Zdata$sampled)

          # random permute to add randomness
          Xdata=Xdata[sample.int(nrow(Xdata)),]
          Zdata=Zdata[sample.int(nrow(Zdata)),]

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
  
      X_sel <- is.na(Xdata[,pair_id])
      Z_sel <- is.na(Zdata[,pair_id])
      Xcounts_single = as.numeric(stats::xtabs(~ factor(Xtype, 1:num_Xu), data=Xdata, subset=X_sel)) 
      Zcounts_single = as.numeric(stats::xtabs(~ factor(Ztype, 1:num_Zu), data=Zdata, subset=Z_sel))
  
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
      
      # The number of people in the population
      N = num_women + num_men
      gw = log(num_women/N) # to ensure exp(gw)+exp(gm) = 1
      gm = log(num_men/N) # to ensure exp(gw)+exp(gm) = 1

      pmf = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(pmf) <- c(cnM,"singles")
      rownames(pmf) <- c(cnW,"singles")
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
  
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
        counts[1:num_Xu,1:num_Zu] <- 2*counts[1:num_Xu,1:num_Zu]
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
      solution[is.na(solution) | solution>=UB] = UB[is.na(solution) | solution>=UB]-0.01
      solution[is.na(solution) | solution<=LB] = LB[is.na(solution) | solution<=LB]+0.01
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
        message(sprintf("Optimization for starting value %d is overly constrained. Estimates may be unstable.",i))
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


      PMF_SW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
      PMF_SW <- PMF_SW / (PMF_SW + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],1,sum))
      PMF_SM <- pmf_est[ nrow(pmf_est),-ncol(pmf_est),drop=FALSE]
      PMF_SM <- PMF_SM / (PMF_SM + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],2,sum))
      PMF_PW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
      PMF_PW <- 1 - PMF_SW
      PMF_PM <- 1 - PMF_SM 
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

      pmf_est[-nrow(pmf_est), -ncol(pmf_est)] <- 2*pmf_est[-nrow(pmf_est), -ncol(pmf_est)]
      pmf_est <- pmf_est/sum(pmf_est)

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
