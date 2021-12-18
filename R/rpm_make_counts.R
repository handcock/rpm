rpm_make_counts <- function(Xdata, Zdata, sampling_design, sampled, Xid, Zid, pair_id, X_w, Z_w, Xu, Zu, verbose=FALSE){

    num_Xu <- nrow(Xu)
    num_Zu <- nrow(Zu)
    cnW <- paste(colnames(Xu)[2],Xu[,2], sep=".")
    for(i in 2:ncol(Xu)){
      cnW <- paste(cnW,paste(colnames(Xu)[i],Xu[,i], sep="."),sep='.')
    }
    cnM <- paste(colnames(Zu)[2],Zu[,2], sep=".")
    for(i in 2:ncol(Zu)){
      cnM <- paste(cnM,paste(colnames(Zu)[i],Zu[,i], sep="."),sep='.')
    }

    if(sampling_design == "stock-stock"){
    # IDs of the women matched to the sampled men (and vice versa)
     paired_and_sampled_M <- Zdata[,sampled] & !is.na(Zdata[,pair_id])
     paired_and_sampled_W <- Xdata[,sampled] & !is.na(Xdata[,pair_id])
     paired_and_unsampled_M <- !Zdata[,sampled] & !is.na(Zdata[,pair_id])
     paired_and_unsampled_W <- !Xdata[,sampled] & !is.na(Xdata[,pair_id])
     W_paired_to_sampled_M <- match(Zdata[paired_and_sampled_M,pair_id], Xdata[,Xid])
     M_paired_to_sampled_W <- match(Xdata[paired_and_sampled_W,pair_id], Zdata[,Zid])
     W_paired_to_unsampled_M <- match(Zdata[paired_and_unsampled_M,pair_id], Xdata[,Xid])
     M_paired_to_unsampled_W <- match(Xdata[paired_and_unsampled_W,pair_id], Zdata[,Zid])
    }
    if(sampling_design == "stock-flow"){
    # IDs of the women matched to the sampled men (and vice versa)
    # Same as stock-stock code
     paired_and_sampled_M <- Zdata[,sampled] & !is.na(Zdata[,pair_id])
     paired_and_sampled_W <- Xdata[,sampled] & !is.na(Xdata[,pair_id])
     paired_and_unsampled_M <- !Zdata[,sampled] & !is.na(Zdata[,pair_id])
     paired_and_unsampled_W <- !Xdata[,sampled] & !is.na(Xdata[,pair_id])
     W_paired_to_sampled_M <- match(Zdata[paired_and_sampled_M,pair_id], Xdata[,Xid])
     M_paired_to_sampled_W <- match(Xdata[paired_and_sampled_W,pair_id], Zdata[,Zid])
     W_paired_to_unsampled_M <- match(Zdata[paired_and_unsampled_M,pair_id], Xdata[,Xid])
     M_paired_to_unsampled_W <- match(Xdata[paired_and_unsampled_W,pair_id], Zdata[,Zid])
    }
    if(sampling_design == "census"){
     paired_and_sampled_M <- !is.na(Zdata[,pair_id])
     paired_and_sampled_W <- !is.na(Xdata[,pair_id])
     W_paired_to_sampled_M <- match(Zdata[paired_and_sampled_M,pair_id], Xdata[,Xid])
     M_paired_to_sampled_W <- match(Xdata[paired_and_sampled_W,pair_id], Zdata[,Zid])
    }

#   if(is.null(pair_w) & sampling_design != "census"){
#   # This is the unweighted case
#   # Construct individual weight case 
#     a=rep(NA, length=nrow(Xdata))
#     a[paired_and_sampled_W] <- Xdata[paired_and_sampled_W,X_w]
#     a[paired_and_unsampled_W] <- 0
#     Xdata[["pair_w"]]=a
#     a=rep(NA, length=nrow(Zdata))
#     a[paired_and_sampled_M] <- Zdata[paired_and_sampled_M,Z_w]
#     a[paired_and_unsampled_M] <- 0
#     Zdata[["pair_w"]]=a
#     pair_w <- "pair_w"
#   }
#   if(is.null(pair_w)){ pair_w <- "pair_w" }
    
    # get the proportion of men and women
    
    if(sampling_design != "census"){
     if (sampling_design == "stock-stock") {
      XdataM <- Xdata[W_paired_to_sampled_M,]
      ZdataW <- Zdata[M_paired_to_sampled_W,]
      # Presumes individual weights
      num_women = sum(Xdata[Xdata[,sampled] & is.na(Xdata[,pair_id]),X_w])
      num_women = num_women + sum(Xdata[Xdata[,sampled] & !is.na(Xdata[,pair_id]),X_w]) # number of women
      num_men = sum(Zdata[Zdata[,sampled] & is.na(Zdata[,pair_id]),Z_w])
      num_men = num_men + sum(Zdata[Zdata[,sampled] & !is.na(Zdata[,pair_id]),Z_w]) # number of men
      n = num_women + num_men
      gw = log(num_women/n) # to ensure exp(gw)+exp(gm) = 1
      gm = log(num_men/n) # to ensure exp(gw)+exp(gm) = 1
     }else{
      # for now this is stock-flow
      # Presumes individual weights
      num_women = sum(Xdata[Xdata[,sampled],X_w])
      num_men = sum(Zdata[Zdata[,sampled],Z_w])
      n = num_women + num_men
      gw = log(num_women/n) # to ensure exp(gw)+exp(gm) = 1
      gm = log(num_men/n) # to ensure exp(gw)+exp(gm) = 1
     }
    }else{
      n = nrow(Xdata) + nrow(Zdata) # The population size
      num_women = nrow(Xdata)
      num_men = nrow(Zdata)
      gw = log(nrow(Xdata)/n) # to ensure exp(gw)+exp(gm) = 1
      gm = log(nrow(Zdata)/n)
    }

    if (sampling_design == "stock-stock") {
     subset=Xdata[,sampled] &  is.na(Xdata[,pair_id])
     pmfW_S = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=Xdata, subset=subset))
     subset=Xdata[,sampled] & !is.na(Xdata[,pair_id])
     pmfW_P = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=Xdata, subset=subset))
     pmfW = pmfW_S + 0.5*pmfW_P
     subset=Zdata[,sampled] &  is.na(Zdata[,pair_id])
     pmfM_S = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=Zdata, subset=subset))
     subset=Zdata[,sampled] & !is.na(Zdata[,pair_id])
     pmfM_P = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=Zdata, subset=subset))
     pmfM = pmfM_S + 0.5*pmfM_P
    }
    if (sampling_design == "stock-flow") {
# What should this be?
     pmfW = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=Xdata, subset=sampled))
     pmfM = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=Zdata, subset=sampled))
    }
    if (sampling_design == "census") {
     pmfW = as.numeric(stats::xtabs(~ factor(Xtype,1:num_Xu), data=Xdata))
     pmfM = as.numeric(stats::xtabs(~ factor(Ztype,1:num_Zu), data=Zdata))
    }
    pmfW = pmfW/sum(pmfW)
    pmfM = pmfM/sum(pmfM)
    names(pmfW) <- cnW
    names(pmfM) <- cnM

    if(verbose){
     message(sprintf("Proportion population paired size: %f",sum(Zdata[paired_and_sampled_M,Z_w])/n))
    }

    if (sampling_design != "census") {
       # These indicate those sampled who are single
       # The counts of the number of women in the population who are single
       #     Xtype_single = as.numeric(stats::xtabs(X_w ~ factor(Xtype, 1:num_Xu), data=Xdata, subset=X_sel))  # account for missing types
       XdataS <- Xdata[Xdata[,sampled],]
       ZdataS <- Zdata[Zdata[,sampled],]
       XdataM <- Xdata[W_paired_to_sampled_M,]
       ZdataW <- Zdata[M_paired_to_sampled_W,]
       X_sel <- is.na(XdataS[,pair_id])
       Z_sel <- is.na(ZdataS[,pair_id])
  
       Xtype_single = as.numeric(stats::xtabs(X_w ~ factor(Xtype, 1:num_Xu), data=XdataS, subset=X_sel)) 
       Ztype_single = as.numeric(stats::xtabs(Z_w ~ factor(Ztype, 1:num_Zu), data=ZdataS, subset=Z_sel))
     } else {
       XdataS <- Xdata
       ZdataS <- Zdata
       XdataM <- Xdata[W_paired_to_sampled_M,]
       ZdataW <- Zdata[M_paired_to_sampled_W,]
       X_sel <- is.na(XdataS[,pair_id])
       Z_sel <- is.na(ZdataS[,pair_id])
    }
    # Compute the number sampled (mainly for the s.e. computation)
    if (sampling_design == "stock-stock") {
       num_sampled <- nrow(Xdata)+nrow(Zdata)-sum(!is.na(Zdata[,pair_id]))
    }else{
      if (sampling_design == "stock-flow") {
        num_sampled <- sum(Xdata[,sampled])+sum(Zdata[,sampled])
      }else{
        # census for now
        num_sampled <- nrow(Xdata)+nrow(Zdata)-sum(!is.na(Zdata[,pair_id]))
      }
    }
    Xcounts_single = as.numeric(stats::xtabs(~ factor(Xtype, 1:num_Xu), data=XdataS, subset=X_sel)) 
    Zcounts_single = as.numeric(stats::xtabs(~ factor(Ztype, 1:num_Zu), data=ZdataS, subset=Z_sel))

    if (sampling_design == "census") { 
      
      # The number of people in the population
      N <- n

      pmf = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(pmf) <- c(cnM,"singles")
      rownames(pmf) <- c(cnW,"singles")
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
      
      # The number of people in the population of each (X,Z) pair
      pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      # The proportion of people in the population of each (X,Z) pair
      pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N
      # The proportion of households in the population of each (X,Z) pair

      # The number of pairs of people in the sample of each (X,Z) pair
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~ Xtype_paired + Ztype_paired))
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~XdataS[paired_and_sampled_W,"Xtype"]+ZdataW[,"Ztype"])) + as.numeric(stats::xtabs(~XdataM[,"Xtype"]+ZdataS[paired_and_sampled_M,"Ztype"]))
     #counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
      counts[1:num_Xu,1:num_Zu] <- 2*counts[1:num_Xu,1:num_Zu]
      
      if (!is.empty(Xcounts_single)) {
        pmf[1:num_Xu,1+num_Zu] = Xcounts_single / N
        counts[1:num_Xu,1+num_Zu] = Xcounts_single
      }
      if (!is.empty(Zcounts_single)) {
        pmf[1+num_Xu,1:num_Zu] = Zcounts_single / N
        counts[1+num_Xu,1:num_Zu] = Zcounts_single
      }

      # The proportion of households in the population of each (X,Z) pair
      #pmf <- pmf / sum(pmf)

      pmfN <- pmf*num_sampled
  #   pmfN[1:num_Xu,1:num_Zu] <- pmfN[1:num_Xu,1:num_Zu] / 2
  #   pmfN <- counts 

    } else if (sampling_design == "stock-stock") {
      
      pmf = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(pmf) <- c(cnM,"singles")
      rownames(pmf) <- c(cnW,"singles")
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
      
      # The number of people in the population
      N <- n

      # The number of people in the population of each (X,Z) pair
#     pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XZ_paired_w ~ Xtype_paired + Ztype_paired))
#     pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XdataS[,X_w]~XdataS[,"Xtype"]+ZdataW[,"Ztype"])) + as.numeric(stats::xtabs(ZdataS[,Z_w]~XdataM[,"Xtype"]+ZdataS[,"Ztype"]))
      pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      # The proportion of households in the population of each (X,Z) pair
      #pmf[1:num_Xu,1:num_Zu] <- 0.5*pmf[1:num_Xu,1:num_Zu] / N
      # The proportion of people in the population of each (X,Z) pair
      pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N

      # The number of pairs of people in the sample of each (X,Z) pair
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~ Xtype_paired + Ztype_paired))
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~XdataS[,"Xtype"]+ZdataW[,"Ztype"])) + as.numeric(stats::xtabs(~XdataM[,"Xtype"]+ZdataS[,"Ztype"]))
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
#     counts[1:num_Xu,1:num_Zu] <- 0.5*counts[1:num_Xu,1:num_Zu]
      counts[1:num_Xu,1:num_Zu] <- 2*counts[1:num_Xu,1:num_Zu]
      
      if (!is.empty(Xtype_single)) {
        pmf[1:num_Xu,1+num_Zu] = Xtype_single / N
        counts[1:num_Xu,1+num_Zu] = Xcounts_single
      }
      if (!is.empty(Ztype_single)) {
        pmf[1+num_Xu,1:num_Zu] = Ztype_single / N
        counts[1+num_Xu,1:num_Zu] = Zcounts_single
      }

      # The proportion of households in the population of each (X,Z) pair
      # pmf <- pmf / sum(pmf)

#     Next needed for microsimulate
     #pmfN[1:num_Xu,1:num_Zu] <- pmfN[1:num_Xu,1:num_Zu] / 2
#     Next needed for simulate
      pmfN <- counts
    # pmfN[1:num_Xu,1:num_Zu] <- pmfN[1:num_Xu,1:num_Zu] / 2
    # pmfN <- pmf*num_sampled

    } else { # assume "stock-flow"
      
      pmf = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(pmf) <- c(cnM,"singles")
      rownames(pmf) <- c(cnW,"singles")
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
      #pmf_alt <- pmf
      
      # The number of people in the population
      N <- n

      # The number of people in the population of each (X,Z) pair
#     pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XZ_paired_w ~ Xtype_paired + Ztype_paired))
#     pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XdataS[,X_w]~XdataS[,"Xtype"]+ZdataW[,"Ztype"])) + as.numeric(stats::xtabs(ZdataS[,Z_w]~XdataM[,"Xtype"]+ZdataS[,"Ztype"]))
# orig      pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(Zdata[M_paired_to_sampled_W,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(Xdata[W_paired_to_sampled_M,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N
      #pmf_alt[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(ZdataW[,Z_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(XdataM[,X_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      #pmf_alt[1:num_Xu,1:num_Zu] <- 0.5*pmf_alt[1:num_Xu,1:num_Zu] / N
      #pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] + pmf_alt[1:num_Xu,1:num_Zu]
      
      # The proportion of people in the population of each (X,Z) pair
      # note that pmf[1:num_Xu,1:num_Zu] is *twice* the number of pairs so that
      # exp(-gw)*apply(pmf,1,sum) = pmfW
      # pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N
      # The proportion of households in the population of each (X,Z) pair
      # pmf[1:num_Xu,1:num_Zu] <- 0.5*pmf[1:num_Xu,1:num_Zu] / N

      # The number of pairs of people in the sample of each (X,Z) pair
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~ Xtype_paired + Ztype_paired))
#     counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~XdataS[,"Xtype"]+ZdataW[,"Ztype"])) + as.numeric(stats::xtabs(~XdataM[,"Xtype"]+ZdataS[,"Ztype"]))
#
#     Seems twice the stock-stock case
#     so try the same as stock-stock
# orig  counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(Zdata[M_paired_to_sampled_W,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(Xdata[W_paired_to_sampled_M,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))

    # counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
#     counts[1:num_Xu,1:num_Zu] <- 0.5*counts[1:num_Xu,1:num_Zu]
      counts[1:num_Xu,1:num_Zu] <- 2*counts[1:num_Xu,1:num_Zu]
      
      # # account for types of pairs that are not observed
      # # missing rows
      # idx = which(!((1:num_Xu) %in% rownames(tmp)))
      # if (!is.empty(idx)) {
      #   for (ii in 1:length(idx)) {
      #     pmf[(idx[ii]+1):(nrow(tmp)+ii),1:ncol(tmp)]= pmf[idx[ii]:(nrow(tmp)+ii-1),1:ncol(tmp)]
      #     pmf[idx[ii],] = 0
      #   }
      # }
      # # missing cols
      # idx = which(!((1:num_Zu) %in% colnames(tmp)))
      # if (!is.empty(idx)) {
      #   for (ii in 1:length(idx)) {
      #     pmf[,(idx[ii]+1):(ncol(tmp)+ii)]= pmf[,idx[ii]:(ncol(tmp)+ii-1)]
      #     pmf[,idx[ii]] = 0
      #   }
      # }

      if (!is.empty(Xtype_single)) {
        pmf[1:num_Xu,1+num_Zu] = Xtype_single / N
        counts[1:num_Xu,1+num_Zu] = Xcounts_single
      }
      if (!is.empty(Ztype_single)) {
        pmf[1+num_Xu,1:num_Zu] = Ztype_single / N
        counts[1+num_Xu,1:num_Zu] = Zcounts_single
      }
      
      if(verbose){
        message(sprintf("Population size: %f",N))
        message(sprintf("Matrix of counts:"))
        print(counts)
      }

      # The proportion of households in the population of each (X,Z) pair
      # pmf <- pmf / sum(pmf)

#     Next needed for microsimulate
      #pmfN[1:num_Xu,1:num_Zu] <- pmfN[1:num_Xu,1:num_Zu] / 2
#     Next needed for simulate
      pmfN <- counts
    # pmfN[1:num_Xu,1:num_Zu] <- pmfN[1:num_Xu,1:num_Zu] / 2
    # pmfN <- pmf*num_sampled

    }
    if(verbose){
      message(sprintf("Population proportions of women by category:"))
      print(pmfW)
      message(sprintf("Population proportions of men by category:"))
      print(pmfM)
      message(sprintf("Matrix of population households proportions of women x men by category:"))
      print(pmf)
    }
   
    list(pmf=pmf, counts=counts, pmfW=pmfW, pmfM=pmfM, pmfN=pmfN, N=N, gw=gw, gm=gm,
         num_women=num_women, num_men=num_men, num_sampled=num_sampled
        )
}
