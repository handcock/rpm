rpm.bootstrap.small <- function(i, solution, num_women, num_men, Jw, Jm, U_star, V_star,
     S, X, Z, pmfW, pmfM,
     Xu, Zu, num_Xu, num_Zu, cnW, cnM, Xid, Zid, pair_id, sampled, sampling_design,
     NumBeta,NumGammaW,NumGammaM,LB,UB,control,
     num_sampled){

        Ws <- rep(seq_along(pmfW),round(pmfW*num_women))
        if(length(Ws) > num_women+0.01){
         Ws <- Ws[-sample.int(length(Ws),size=length(Ws)-num_women,replace=FALSE)]
        }
        if(length(Ws) < num_women-0.01){
         Ws <- c(Ws,sample.int(length(pmfW),size=num_women-length(Ws),prob=pmfW,replace=TRUE))
        }
        Ws <- Ws[sample.int(length(Ws))]
        Ms <- rep(seq_along(pmfM),round(pmfM*num_men))
        if(length(Ms) > num_men+0.01){
         Ms <- Ms[-sample.int(length(Ms),size=length(Ms)-num_men,replace=FALSE)]
        }
        if(length(Ms) < num_men-0.01){
         Ms <- c(Ms,sample.int(length(pmfM),size=num_men-length(Ms),prob=pmfM,replace=TRUE))
        }
        Ms <- Ms[sample.int(length(Ms))]

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
          # Do not run for a census
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

      mc <- rpm_make_counts(Xdata, Zdata, sampling_design, sampled, Xid, Zid, pair_id, X_w, Z_w, Xu, Zu, verbose=FALSE)
      pmf <- mc$pmf
      num_sampled=mc$num_sampled
      nobs <- num_sampled
      counts=mc$counts
      pmfW=mc$pmfW
      pmfM=mc$pmfM
      pmfN=mc$pmfN
      N=mc$N
      gw=mc$gw
      gm=mc$gm
      num_women=mc$num_women
      num_men=mc$num_men

      control$xtol_rel=control$bs.xtol_rel
      control$maxeval=control$bs.maxeval
      solution[is.na(solution) | solution>=UB] = UB[is.na(solution) | solution>=UB]-0.01
      solution[is.na(solution) | solution<=LB] = LB[is.na(solution) | solution<=LB]+0.01
      out.text <- capture.output(
       out.fit <- nloptr::nloptr(x0=solution, eval_f=loglikfun, 
                 eval_grad_f=gloglikfun,
                 eval_g_eq=eqfun,  eval_jac_g_eq=jeqfun,
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
      pmf_est <- exp(logpmfest(th_hat[1:NumBeta],
                  GammaW=th_hat[NumBeta+(1:NumGammaW)], 
                  GammaM=th_hat[(NumBeta+NumGammaW)+(1:NumGammaM)],
                  S, X, Z,
                  pmfW, pmfM, gw=gw, gm=gm))
      pmf_est[nrow(pmf_est),ncol(pmf_est)] <- 0


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

      if(!control$bs.save.data){Xdata <- Zdata <- NULL}
      list(est=th_hat, LOGODDS_SW=LOGODDS_SW,LOGODDS_SM=LOGODDS_SM,
           Xdata=Xdata,
           Zdata=Zdata,
           pmfW=pmfW, pmfM=pmfM,
           pmf=pmf, counts=counts, nobs=num_sampled,
           pmf_est=pmf_est
          )
     }
