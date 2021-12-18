#' Calculate goodness-of-fit statistics for Revealed Preference Matchings Model based on observed data
#' 
#' \code{\link{gof.rpm}} ...
#' It is typically based on the estimate from a \code{rpm()} call.
#' 
#' The function \code{\link{rpm}} is used to fit a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using an approximate likelihood based on ideas from Menzel (2015).
#' 
#' The model represents the dyadic utility functions as deterministic linear utility functions of
#' dyadic variables. These utility functions are functions of observed characteristics of the women
#' and men.
#' These functions are entered as terms in the function call
#' to \code{\link{rpm}}. This function simulates from such a model.
#'
#' @aliases gof.default
#' @param object list; an object of class\code{rpm} that is typically the result of a call to \code{rpm()}.
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' @param empirical_p logical; (Optional) If TRUE the function returns the empirical p-value of the sample
#' statistic based on \code{nsim} simulations
#' @param nsim integer; (Optional) Number of samples to draw from a large population, simulated based on the model specified
#' by \code{object}; Must be specified if \code{empirical_p} is TRUE
#' @param seed integer; (Optional) random number seed.
#' @param parallel logical, with default FALSE; if `TRUE`, then the function will use parallel processing
#' when computing empirical p value
#' @param ncores integer; number of cores to use for parallel processing; must be specified if
#' parallel is TRUE
#' @param parallel.type string; type of clister to create. Typically, "POCK" or "MPI".
#' parallel is TRUE
#' @param compare_sim string; describes which two objects are compared to compute simulated goodness-of-fit
#' statistics; valid values are \code{"sim-est"}: compares the marginal distribution of pairings in a
#' simulated sample to the \code{rpm} model estimate of the marginal distribution based on that same simulated sample;
#' \code{mod-est}: compares the marginal distribution of pairings in a
#' simulated sample to the \code{rpm} model estimate used to generate the sample
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including data summary statistics.
#' @return \code{\link{gof.rpm}} returns a list consisting of the following elements:
#' \item{obs_pmf}{numeric matrix giving observed probability mass distribution over different household types}
#' \item{observed_pmf_est}{numeric matrix giving expected probability mass distribution from \code{rpm} model}
#' \item{obs_chi_sq}{the count-based observed chi-square statistic comparing marginal distributions of the population
#' the data and the model estimate}
#' \item{obs_chi_sq_cell}{the contribution to the observed chi-squared statistic by household type}
#' \item{obs_kl}{the Kullback-Leibler (KL) divergence computed by comparing the observed marginal distributions to the
#' expected marginal distribution based on the \code{rpm} model estimate}
#' \item{obs_kl_cell}{the contribution to the observed KL divergence by household type}
#' \item{empirical_p_chi_sq}{the proportion of simulated chi-square statistics that are greater than
#' or equal to the observed chi-square statistic}
#' \item{empirical_p_kl}{the proportion of simulated KL divergences that are greater than
#' or equal to the observed KL divergence}
#' \item{chi_sq_simulated}{vector of size \code{nsim} storing all simulated chi-square statistics}
#' \item{kl_simulated}{vector of size \code{nsim} storing all simulated KL divergences}
#' \item{chi_sq_cell_mean}{Mean contributions of each household type to the simulated chi_sq statistic}
#' \item{chi_sq_cell_sd}{Standard deviation of the contributions of each household type to the simulated chi_sq statistics}
#' \item{chi_sq_cell_median}{Median contributions of each household type to the simulated chi_sq statistic}
#' \item{chi_sq_cell_iqr}{Interquartile range of the contributions of each household type to the simulated chi_sq statistics}
#' \item{kl_cell_mean}{Mean contributions of each household type to the simulated KL divergences}
#' \item{kl_cell_sd}{Standard deviation of the contributions of each household type to the simulated KL divergencesc}
#' \item{kl_cell_median}{Median contributions of each household type to the simulated KL divergences}
#' \item{kl_cell_iqr}{Interquartile range of the contributions of each household type to the simulated KL divergences}
#' @keywords models
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' a <- gof(fit)
#' @references Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @export gof
gof <- function(object, ...){
      UseMethod("gof")
    }


#' @noRd
#' @importFrom utils methods
#' @export
gof.default <- function(object,...) {
  classes <- setdiff(gsub(pattern="^gof.",replacement="",as.vector(methods("gof"))), "default")
  stop("Goodness-of-Fit methods have been implemented only for class 'rpm'.")
}

#' @describeIn gof Calculate goodness-of-fit statistics for Revealed Preference Matchings Model based on observed data
#' @export
gof.rpm <- function(object, ...,
  empirical_p = FALSE, nsim=1000, seed = NULL, parallel = FALSE, ncores = NA, parallel.type="SOCK", compare_sim = 'sim-est',verbose=FALSE) 
{
  if(!is.null(seed)) set.seed(seed)
  if (empirical_p & is.null(nsim)) {
    stop("When empirical_p is TRUE, please specify the number of simulations used to calculate the empirical p value.")
  }
  if (empirical_p & parallel & is.null(ncores)) {
    stop("For parallel processing, please specify the number of cores in the cluster before re-running.")
  }
  
  observed_pmf_est <- object$pmf_est
  observed_pmf <- object$pmf

  constraints <- match.arg(object$control[["constraints"]], c("none","M_single"))
  constraints <- match(constraints, c("M_single","none")) - 1
  
  out <- list()
  
  if (object$sampling_design == "stock-stock") {
    num_sampled = sum(object$Xdata[,object$sampled] & is.na(object$Xdata[,object$pair_id]))
    num_sampled = num_sampled + sum(object$Zdata[,object$sampled])
  }else{
    if (object$sampling_design == "stock-flow") {
      num_sampled <- sum(object$Xdata[,object$sampled])+sum(object$Zdata[,object$sampled])
    }else{
      num_sampled <- nrow(object$Xdata) + nrow(object$Zdata)
    }
  }

  out$compare_sim <- match.arg(compare_sim, c("sim-est","mod-est"))
  out$obs_pmf <- object$pmf
  out$observed_pmf_est <- object$pmf_est
  
  # Observed chi-square statistic
  compute_chi_sq <- function(observed, expected) {
    csc <- (observed-expected)^2/(observed)
    csc[observed<1e-10] <- 0
    csc[nrow(observed),ncol(observed)] = 0
    return(csc)
  }
  out$obs_chi_sq_cell <- num_sampled*compute_chi_sq(observed_pmf,observed_pmf_est)
  out$obs_chi_sq <- sum(out$obs_chi_sq_cell)
  
  # Observed KL divergence
  compute_kl <- function(observed, expected) {
    klc <- -observed*log(expected/observed)
    klc[observed<1e-10] <- 0
    klc[nrow(observed),ncol(observed)] = 0
    return(klc)
  }
  out$obs_kl_cell <- num_sampled*compute_kl(observed_pmf,observed_pmf_est)
  out$obs_kl <- sum(out$obs_kl_cell)
  
  if (empirical_p) {
  
    Xdata <- object$Xdata
    Zdata <- object$Zdata

    # Reweight to model
    numX <- nrow(object$pmf)-1
    numZ <- ncol(object$pmf)-1
    XdataS0 <- Xdata[Xdata[,object$sampled],]
    ZdataS0 <- Zdata[Zdata[,object$sampled],]
    X_w_new <- XdataS0[,object$X_w]
    Z_w_new <- ZdataS0[,object$Z_w]
    paired_and_sampled_W <- !is.na(XdataS0[,object$pair_id]) 
    paired_and_sampled_M <- !is.na(ZdataS0[,object$pair_id])
    single_and_sampled_W <-  is.na(XdataS0[,object$pair_id]) 
    single_and_sampled_M <-  is.na(ZdataS0[,object$pair_id])
    # Rescale the weights of $pmf to $pmf_est
    for( j in 1:numX){
      a <- paired_and_sampled_W & XdataS0$Xtype==j
      b <- match(XdataS0[a,object$pair_id], Zdata[,object$Zid])
      w <- XdataS0[a,object$X_w]
      for( k in 1:numZ){
        d <- Zdata$Ztype[b] == k
        if(any(d)){
          w[d] <- w[d] * observed_pmf_est[j,k] / object$pmf[j,k]
        }
      }
      X_w_new[a] <- w
    }
    for( k in 1:numZ){
      a <- paired_and_sampled_M & ZdataS0$Ztype==k
      b <- match(ZdataS0[a,object$pair_id], Xdata[,object$Xid])
      w <- ZdataS0[a,object$Z_w]
      for( j in 1:numX){
        d <- Xdata$Xtype[b] == j
        if(any(d)){
          w[d] <- w[d] * observed_pmf_est[j,k] / object$pmf[j,k]
        }
      }
      Z_w_new[a] <- w
    }
    w <- XdataS0[single_and_sampled_W,object$X_w]
    for( j in 1:numX){
      a <- XdataS0$Xtype[single_and_sampled_W]==j
      if(any(a)){
        w[a] <- w[a] * observed_pmf_est[j,ncol(object$pmf)] / object$pmf[j,ncol(object$pmf)]
      }
    }
    X_w_new[single_and_sampled_W] <- w
    
    w <- ZdataS0[single_and_sampled_M,object$Z_w]
    for( k in 1:numZ){
      a <- ZdataS0$Ztype[single_and_sampled_M]==k
      if(any(a)){
        w[a] <- w[a] * observed_pmf_est[nrow(object$pmf),k] / object$pmf[nrow(object$pmf),k]
      }
    }
    Z_w_new[single_and_sampled_M] <- w

    Xdata[Xdata[,object$sampled],object$X_w] <- X_w_new
    Zdata[Zdata[,object$sampled],object$Z_w] <- Z_w_new
    
    X_w_rel <- rep(1,nrow(XdataS0))
    Z_w_rel <- rep(1,nrow(ZdataS0))
    if(object$sampling_design == "stock-stock"){
      paired_W <- !is.na(XdataS0[,object$pair_id])
      paired_M <- !is.na(ZdataS0[,object$pair_id])
      X_w_rel[paired_W] <- 1.0
      Z_w_rel[paired_M] <- 0.0
    }

    Xu <- unique(Xdata$Xtype)
    Xu <- Xu[do.call(order, as.data.frame(Xu))]
    Zu <- unique(Zdata$Ztype)
    Zu <- Zu[do.call(order, as.data.frame(Zu))]
    
    num_Xu = length(Xu)
    num_Zu = length(Zu)
    
    nbootstrap = nsim

    if (parallel) {
      doFuture::registerDoFuture()
      cl <- parallel::makeCluster(object$control$ncores)
      future::plan(cluster, workers = cl)
      if(Sys.info()[["sysname"]] == "Windows"){
        future::plan(multisession)  ## on MS Windows
      }else{
        future::plan(multicore)     ## on Linux, Solaris, and macOS
      }
      if (!is.null(object$control$seed)) {
        doRNG::registerDoRNG(object$control$seed)
      }
      
      pmf_sim_all <-
        foreach::foreach (b=1:nbootstrap, .packages=c('rpm')
        ) %dorng% {
          sampling_design = object$sampling_design
          sampled = object$sampled
          
          pmfboot = matrix(0, nrow = nrow(observed_pmf), ncol = ncol(observed_pmf))
          colnames(pmfboot) <- colnames(observed_pmf)
          rownames(pmfboot) <- rownames(observed_pmf)
          countsboot <- pmfboot
          
          I <- sample(c(XdataS0[,object$Xid], ZdataS0[,object$Zid]),
                      replace=TRUE,size=num_sampled,
                      prob=c(X_w_rel, Z_w_rel))
#         I <- sample(I, replace=TRUE,size=num_sampled) # Finite population adjustment

          Xmatch <- match(I,Xdata[,object$Xid])
          Zmatch <- match(I,Zdata[,object$Zid])
          XdataS <- Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataS <- Zdata[Zmatch[!is.na(Zmatch)],]

#         Set the sampled indicator
          XdataS[,object$sampled] <- TRUE
          ZdataS[,object$sampled] <- TRUE
          
          # Find the people paired to the sampled people
          paired_and_sampled_W <- !is.na(XdataS[,object$pair_id])
          paired_and_sampled_M <- !is.na(ZdataS[,object$pair_id])
          W_paired_to_sampled_M <- match(ZdataS[paired_and_sampled_M,object$pair_id], Xdata[,object$Xid])
          M_paired_to_sampled_W <- match(XdataS[paired_and_sampled_W,object$pair_id], Zdata[,object$Zid])
          XdataM <- Xdata[W_paired_to_sampled_M,]
          ZdataW <- Zdata[M_paired_to_sampled_W,]
          
          if(object$sampling_design == "stock-flow"){
            XdataM[,object$sampled] <- FALSE
            ZdataW[,object$sampled] <- FALSE
          }else{
        #   XdataM[,object$sampled] <- TRUE
            ZdataW[,object$sampled] <- TRUE
          }

#         #  Merge the sampled and paired-to-sampled
#         XdataB <- rbind(XdataS, XdataM)
#         ZdataB <- rbind(ZdataS, ZdataW)

          X_sel <- is.na(XdataS[,object$pair_id])
          Z_sel <- is.na(ZdataS[,object$pair_id])
          
          Xcounts_single = as.numeric(stats::xtabs(~ factor(Xtype, 1:num_Xu), data=XdataS, subset=X_sel)) 
          Zcounts_single = as.numeric(stats::xtabs(~ factor(Ztype, 1:num_Zu), data=ZdataS, subset=Z_sel))
          Xtype_single = as.numeric(stats::xtabs(X_w ~ factor(Xtype, 1:num_Xu), data=XdataS, subset=X_sel))
          Ztype_single = as.numeric(stats::xtabs(Z_w ~ factor(Ztype, 1:num_Zu), data=ZdataS, subset=Z_sel))
          
          # The number of people in the population
          if (sampling_design == "census") {
            pmfW = as.numeric(stats::xtabs(~ factor(Xtype,1:num_Xu), data=XdataS))
            pmfM = as.numeric(stats::xtabs(~ factor(Ztype,1:num_Zu), data=ZdataS))
            N = nrow(XdataS) + nrow(ZdataS) # The population size
            gw = log(nrow(XdataS)/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(nrow(ZdataS)/N)
          }
          if (sampling_design == "stock-stock") {
            N_w = sum(XdataS[is.na(XdataS[,object$pair_id]),object$X_w])
            N_w = N_w + sum(ZdataS[!is.na(ZdataS[,object$pair_id]),object$Z_w]) # The population size
            N_m = sum(ZdataS[is.na(ZdataS[,object$pair_id]),object$Z_w])
            N_m = N_m + sum(XdataS[!is.na(XdataS[,object$pair_id]),object$X_w]) # The population size
            N = N_w + N_m
            gw = log(N_w/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(N_m/N) # to ensure exp(gw)+exp(gm) = 1
            #
            subset=is.na(XdataS[,object$pair_id])
            pmfW_S = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=XdataS, subset=subset))
            subset=!is.na(XdataS[,object$pair_id])
            pmfW_P = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=XdataS, subset=subset))
            pmfW = pmfW_S + 0.5*pmfW_P
            subset=is.na(ZdataS[,object$pair_id])
            pmfM_S  = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=ZdataS, subset=subset))
            subset=!is.na(ZdataS[,object$pair_id])
            pmfM_P = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=ZdataS, subset=subset))
            pmfM = pmfM_S + 0.5*pmfM_P
          }
          if (sampling_design == "stock-flow") {
            N = sum(XdataS[,object$X_w]) + sum(ZdataS[,object$Z_w]) # The population size
            gw = log(sum(XdataS[,object$X_w])/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(sum(ZdataS[,object$Z_w])/N)

            pmfW = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=XdataS))
            pmfM = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=ZdataS))
          }
          pmfW = pmfW/sum(pmfW)
          pmfM = pmfM/sum(pmfM)

          if (sampling_design == "census") {
            pmfboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
            countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
          }
          if (sampling_design == "stock-stock") {
            pmfboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XdataS[paired_and_sampled_W,object$X_w]~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(ZdataS[paired_and_sampled_M,object$Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
            countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
          }
          if (sampling_design == "stock-flow") {
            pmfboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XdataS[paired_and_sampled_W,object$X_w]~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(ZdataS[paired_and_sampled_M,object$Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
          # countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
            countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
          }

          pmfboot[1:num_Xu,1:num_Zu] <- pmfboot[1:num_Xu,1:num_Zu] / N

          if (!is.empty(Xtype_single)) {
            pmfboot[1:num_Xu,1+num_Zu] = Xtype_single / N
            countsboot[1:num_Xu,1+num_Zu] = Xcounts_single
          }
          if (!is.empty(Ztype_single)) {
            pmfboot[1+num_Xu,1:num_Zu] = Ztype_single / N
            countsboot[1+num_Xu,1:num_Zu] = Zcounts_single
          }
          
          pmfbootN <- pmfboot*num_sampled

          bout <- nloptr::nloptr(x0=object$coefficients, eval_f=loglikfun_default,
                 eval_grad_f=gloglikfun_default,
                 eval_g_eq=eqfun_default, eval_jac_g_eq=jeqfun_default,
                 lb=c(rep(-6,object$NumBeta),rep(-6,object$NumGamma)),
                 ub=c(rep( 6,object$NumBeta),rep( 6,object$NumGamma)),
                 Sd=object$Sd,Xd=object$Xd,Zd=object$Zd,
                 NumGammaW=object$NumGammaW, NumGammaM=object$NumGammaM,
                 pmfW=pmfW, pmfM=pmfM,
                 pmf=pmfboot, counts=pmfbootN, gw=gw, gm=gm, N=N,
                 sampling=object$sampling_design, constraints=constraints,
                 opts=object$control)

           pmfboot_est <- exp(augpmf(bout$coefficients[1:object$NumBeta],
                    GammaW=bout$coefficients[object$NumBeta+(1:object$NumGammaW)],
                    GammaM=bout$coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
                    Sd=object$Sd,Xd=object$Xd,Zd=object$Zd,
                    Sdim=dim(object$Sd), Xdim=dim(object$Xd),Zdim=dim(object$Zd),
                    pmfW=pmfW, pmfM=pmfM,
                    gw=gw, gm=gm))
#          if (sampling_design == "stock-stock") {
#            a <- sum(pmfboot_est[-nrow(pmfboot_est),-ncol(pmfboot_est)])
#            pmfboot_est[-nrow(pmfboot_est),-ncol(pmfboot_est)] <- 2*pmfboot_est[-nrow(pmfboot_est),-ncol(pmfboot_est)]
#            pmfboot_est[-nrow(pmfboot_est), ncol(pmfboot_est)] <- pmfboot_est[-nrow(pmfboot_est), ncol(pmfboot_est)]*(1-2*a)/(1-a)
#            pmfboot_est[ nrow(pmfboot_est),-ncol(pmfboot_est)] <- pmfboot_est[ nrow(pmfboot_est),-ncol(pmfboot_est)]*(1-2*a)/(1-a)
#          }
           pmfboot_est[nrow(pmfboot_est),ncol(pmfboot_est)] <- 0

          return(list(countsboot, pmfboot, pmfboot_est, N, num_sampled))
        }
        parallel::stopCluster(cl)
      } else { # not parallel
        pmf_sim_all <- list()
        for (it in 1:nbootstrap) {

          sampling_design = object$sampling_design
          sampled = object$sampled
          
          pmfboot = matrix(0, nrow = nrow(observed_pmf), ncol = ncol(observed_pmf))
          colnames(pmfboot) <- colnames(observed_pmf)
          rownames(pmfboot) <- rownames(observed_pmf)
          countsboot <- pmfboot
          
          I <- sample(c(XdataS0[,object$Xid], ZdataS0[,object$Zid]),
                      replace=TRUE,size=num_sampled,
                      prob=c(X_w_rel, Z_w_rel))
          I <- sample(I, replace=TRUE,size=num_sampled) # Finite population adjustment

          Xmatch <- match(I,Xdata[,object$Xid])
          Zmatch <- match(I,Zdata[,object$Zid])
          XdataS <- Xdata[Xmatch[!is.na(Xmatch)],]
          ZdataS <- Zdata[Zmatch[!is.na(Zmatch)],]
          
          # Find the people paired to the sampled people
          paired_and_sampled_W <- !is.na(XdataS[,object$pair_id])
          paired_and_sampled_M <- !is.na(ZdataS[,object$pair_id])
          W_paired_to_sampled_M <- match(ZdataS[paired_and_sampled_M,object$pair_id], Xdata[,object$Xid])
          M_paired_to_sampled_W <- match(XdataS[paired_and_sampled_W,object$pair_id], Zdata[,object$Zid])
          XdataM <- Xdata[W_paired_to_sampled_M,]
          ZdataW <- Zdata[M_paired_to_sampled_W,]
          
          X_sel <- is.na(XdataS[,object$pair_id])
          Z_sel <- is.na(ZdataS[,object$pair_id])
          
          Xcounts_single = as.numeric(stats::xtabs(~ factor(Xtype, 1:num_Xu), data=XdataS, subset=X_sel)) 
          Zcounts_single = as.numeric(stats::xtabs(~ factor(Ztype, 1:num_Zu), data=ZdataS, subset=Z_sel))
          Xtype_single = as.numeric(stats::xtabs(X_w ~ factor(Xtype, 1:num_Xu), data=XdataS, subset=X_sel))
          Ztype_single = as.numeric(stats::xtabs(Z_w ~ factor(Ztype, 1:num_Zu), data=ZdataS, subset=Z_sel))
          
          # The number of people in the population
          if (sampling_design == "census") {
            pmfW = as.numeric(stats::xtabs(~ factor(Xtype,1:num_Xu), data=XdataS))
            pmfM = as.numeric(stats::xtabs(~ factor(Ztype,1:num_Zu), data=ZdataS))
            N = nrow(XdataS) + nrow(ZdataS) # The population size
            gw = log(nrow(XdataS)/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(nrow(ZdataS)/N)
          }
          if (sampling_design == "stock-stock") {
            N_w = sum(XdataS[is.na(XdataS[,object$pair_id]),object$X_w])
            N_w = N_w + sum(ZdataS[!is.na(ZdataS[,object$pair_id]),object$Z_w]) # The population size
            N_m = sum(ZdataS[is.na(ZdataS[,object$pair_id]),object$Z_w])
            N_m = N_m + sum(XdataS[!is.na(XdataS[,object$pair_id]),object$X_w]) # The population size
            N = N_w + N_m
            gw = log(N_w/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(N_m/N) # to ensure exp(gw)+exp(gm) = 1
            #
            subset=is.na(XdataS[,object$pair_id])
            pmfW_S = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=XdataS, subset=subset))
            subset=!is.na(XdataS[,object$pair_id])
            pmfW_P = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=XdataS, subset=subset))
            pmfW = pmfW_S + 0.5*pmfW_P
            subset=is.na(ZdataS[,object$pair_id])
            pmfM_S  = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=ZdataS, subset=subset))
            subset=!is.na(ZdataS[,object$pair_id])
            pmfM_P = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=ZdataS, subset=subset))
            pmfM = pmfM_S + 0.5*pmfM_P
          }
          if (sampling_design == "stock-flow") {
            N = sum(XdataS[,object$X_w]) + sum(ZdataS[,object$Z_w]) # The population size
            gw = log(sum(XdataS[,object$X_w])/N) # to ensure exp(gw)+exp(gm) = 1
            gm = log(sum(ZdataS[,object$Z_w])/N)

            pmfW = as.numeric(stats::xtabs(X_w ~ factor(Xtype,1:num_Xu), data=XdataS))
            pmfM = as.numeric(stats::xtabs(Z_w ~ factor(Ztype,1:num_Zu), data=ZdataS))
          }
          pmfW = pmfW/sum(pmfW)
          pmfM = pmfM/sum(pmfM)

          if (sampling_design == "census") {
            pmfboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
            countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
          }
          if (sampling_design == "stock-stock") {
            pmfboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XdataS[paired_and_sampled_W,object$X_w]~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(ZdataS[paired_and_sampled_M,object$Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
            countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
          }
          if (sampling_design == "stock-flow") {
            pmfboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(XdataS[paired_and_sampled_W,object$X_w]~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(ZdataS[paired_and_sampled_M,object$Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
        #   countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
            countsboot[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(XdataS[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(ZdataS[paired_and_sampled_M,"Ztype"],1:num_Zu)))
          }

          pmfboot[1:num_Xu,1:num_Zu] <- pmfboot[1:num_Xu,1:num_Zu] / N
          pmfboot[1:num_Xu,1+num_Zu] = Xtype_single / N
          pmfboot[1+num_Xu,1:num_Zu] = Ztype_single / N

          if (!is.empty(Xcounts_single)) {
            countsboot[1:num_Xu,1+num_Zu] = Xcounts_single
          }
          if (!is.empty(Zcounts_single)) {
            countsboot[1+num_Xu,1:num_Zu] = Zcounts_single
          }
          
          pmfbootN <- pmfboot*num_sampled

          bout <- nloptr::nloptr(x0=object$coefficients, eval_f=loglikfun,
                 eval_grad_f=gloglikfun,
                 eval_g_eq=eqfun, eval_jac_g_eq=jeqfun,
                 lb=c(rep(-6,object$NumBeta),rep(-6,object$NumGamma)),
                 ub=c(rep( 6,object$NumBeta),rep( 6,object$NumGamma)),
                 Sd=object$Sd,Xd=object$Xd,Zd=object$Zd,
                 NumGammaW=object$NumGammaW, NumGammaM=object$NumGammaM,
                 pmfW=pmfW, pmfM=pmfM,
                 pmf=pmfboot, counts=pmfbootN, gw=gw, gm=gm, N=N,
                 sampling=object$sampling_design, constraints=constraints,
                 opts=object$control)

           pmfboot_est <- exp(augpmf(bout$coefficients[1:object$NumBeta],
                    GammaW=bout$coefficients[object$NumBeta+(1:object$NumGammaW)],
                    GammaM=bout$coefficients[(object$NumBeta+object$NumGammaW)+(1:object$NumGammaM)],
                    Sd=object$Sd,Xd=object$Xd,Zd=object$Zd,
                    Sdim=dim(object$Sd), Xdim=dim(object$Xd),Zdim=dim(object$Zd),
                    pmfW=pmfW, pmfM=pmfM,
                    gw=gw, gm=gm))
#          if (sampling_design == "stock-stock") {
#            a <- sum(pmfboot_est[-nrow(pmfboot_est),-ncol(pmfboot_est)])
#            pmfboot_est[-nrow(pmfboot_est),-ncol(pmfboot_est)] <- 2*pmfboot_est[-nrow(pmfboot_est),-ncol(pmfboot_est)]
#            pmfboot_est[-nrow(pmfboot_est), ncol(pmfboot_est)] <- pmfboot_est[-nrow(pmfboot_est), ncol(pmfboot_est)]*(1-2*a)/(1-a)
#            pmfboot_est[ nrow(pmfboot_est),-ncol(pmfboot_est)] <- pmfboot_est[ nrow(pmfboot_est),-ncol(pmfboot_est)]*(1-2*a)/(1-a)
#          }
           pmfboot_est[nrow(pmfboot_est),ncol(pmfboot_est)] <- 0
           pmf_sim_all[[it]] <- list(countsboot, pmfboot, pmfboot_est, N, num_sampled)
        }
    }
    
    samp_sizes <- unlist(lapply(pmf_sim_all, '[[', 5))
    simulated_pmf <- array(unlist(lapply(pmf_sim_all, '[[', 2)),dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    simulated_pmf_est <- array(unlist(lapply(pmf_sim_all, '[[', 3)),dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    
    if(verbose){
      message(sprintf("Assessing rpm model goodness-of-fit: %s",out$compare_sim))
    }
    # Calculate chi-square contribution by cell for simulated pmfs
    if (out$compare_sim == 'sim-est') {
      sim_chi_sq_cell <- array(dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      for (it in 1:nsim) {
        sim_chi_sq_cell[,,it] = samp_sizes[it]*compute_chi_sq(simulated_pmf[,,it],simulated_pmf_est[,,it])
      }
      
      # Calculate KL contribution by cell for simulated pmfs
      sim_kl_cell <- array(NA, c(nrow(object$pmf),ncol(object$pmf),nsim))
      for (it in 1:nsim) {
        sim_kl_cell[,,it] <- samp_sizes[it]*compute_kl(simulated_pmf[,,it],simulated_pmf_est[,,it])
      }
    } else {
      sim_chi_sq_cell <- array(dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      for (it in 1:nsim) {
        sim_chi_sq_cell[,,it] = samp_sizes[it]*compute_chi_sq(simulated_pmf[,,it],observed_pmf_est)
      }
      
      # Calculate KL contribution by cell for simulated pmfs
      sim_kl_cell <- array(NA, c(nrow(object$pmf),ncol(object$pmf),nsim))
      for (it in 1:nsim) {
        sim_kl_cell[,,it] <- samp_sizes[it]*compute_kl(simulated_pmf[,,it],observed_pmf_est)
      }
    }
    
    # Simulated chi_sq and KL divergence statistics
    out$chi_sq_simulated <- apply(sim_chi_sq_cell, 3, sum, na.rm=TRUE)
    out$kl_simulated <- apply(sim_kl_cell, 3, sum, na.rm=TRUE)
    
    # Empirical p-value
    out$empirical_p_chi_sq <- mean(out$chi_sq_simulated >= out$obs_chi_sq, na.rm=TRUE)
    out$empirical_p_kl <- mean(out$kl_simulated >= out$obs_kl, na.rm=TRUE)
    
    # Cell summaries of contributions
    # Mean, SD, Median, IQR for chi-sq
    out$chi_sq_cell_mean <- apply(sim_chi_sq_cell, 1:2, mean, na.rm=TRUE)
    dimnames(out$chi_sq_cell_mean) <- dimnames(observed_pmf)

    out$chi_sq_cell_sd <- apply(sim_chi_sq_cell, 1:2, stats::sd, na.rm=TRUE)
    dimnames(out$chi_sq_cell_sd) <- dimnames(observed_pmf)    
    
    out$chi_sq_cell_median <- apply(sim_chi_sq_cell, 1:2, stats::median, na.rm=TRUE)
    dimnames(out$chi_sq_cell_median) <- dimnames(observed_pmf)
    
    out$chi_sq_cell_iqr <- apply(sim_chi_sq_cell, 1:2, stats::IQR, na.rm=TRUE)
    dimnames(out$chi_sq_cell_iqr) <- dimnames(observed_pmf)
    
    # Mean, SD, Median, IQR for KL
    out$kl_cell_mean <- apply(sim_kl_cell, 1:2, mean, na.rm=TRUE)
    dimnames(out$kl_cell_mean) <- dimnames(observed_pmf)
    
    out$kl_cell_sd <- apply(sim_kl_cell, 1:2, stats::sd, na.rm=TRUE)
    dimnames(out$kl_cell_sd) <- dimnames(observed_pmf)
    
    out$kl_cell_median <- apply(sim_kl_cell, 1:2, stats::median, na.rm=TRUE)
    dimnames(out$kl_cell_median) <- dimnames(observed_pmf)
    
    out$kl_cell_iqr <- apply(sim_kl_cell, 1:2, stats::IQR, na.rm=TRUE)
    dimnames(out$kl_cell_iqr) <- dimnames(observed_pmf)
  }
  
  class(out) <- "gofrpm"
  return(out)
}  

#' @method print gofrpm
#' @export
print.gofrpm <- function(x, ...){
  print( mean(x$chi_sq_sim >= x$obs_chi_sq, na.rm=TRUE) )
  print( mean(x$kl_sim >= x$obs_kl, na.rm=TRUE) )
  invisible()
}

###################################################################
# The <plot.gof> function plots the GOF diagnostics for each
# term included in the build of the gof
#
# --PARAMETERS--
#   x          : a gof, as returned by one of the <gof.X>
#                functions
#   ...        : additional par arguments to send to the native R
#                plotting functions
#   cex.axis   : the magnification of the text used in axis notation;
#                default=0.7
#   plotlogodds: whether the summary results should be presented
#                as their logodds; default=FALSE
#   main       : the main title; default="Goodness-of-fit diagnostics"
#   verbose    : this parameter is ignored; default=FALSE
#   normalize.reachibility: whether to normalize the distances in
#                the 'distance' GOF summary; default=FALSE
#
# --RETURNED--
#   NULL
#
###################################################################



#' @describeIn gof
#' 
#' \code{\link{plot.gofrpm}} plots diagnostics such empirical p-value
#' based on chi-square statistics and KL divergences.
#' See \code{\link{rpm}} for more information on these models.
#' 
#' @param x a list, usually an object of class gofrpm
#' @param cex.axis the magnification of the text used in axis notation;
#  default=0.7#'
#' @param main Title for the goodness-of-fit plots.
#' @keywords graphs
#' 
#' @importFrom graphics points lines mtext plot
#' @export

#' @method plot gofrpm
#' @export plot.gofrpm
plot.gofrpm <- function(x, ..., 
         cex.axis=0.7,
         main="Goodness-of-fit diagnostics"
		 ) {

	  hist(x$chi_sq_sim, xlim=range(c(x$chi_sq_sim,x$obs_chi_sq),finite=TRUE),
		 probability=TRUE,
		 main=expression(paste("Distribution of simulated ",chi^2)),
		 xlab=expression(chi^2),nclass=20)
	  lines(stats::density(x$chi_sq_sim),col=3)
	  abline(v=x$obs_chi_sq, col=2)

	  hist(x$kl_sim, xlim=range(c(x$kl_sim,x$obs_kl),finite=TRUE),
		 probability=TRUE,
		 main=expression(paste("Distribution of simulated ",KL)),
		 xlab=expression(KL),nclass=20)
	  lines(stats::density(x$kl_sim),col=3)
	  abline(v=x$obs_kl, col=2)
	  
	y <- x$chi_sq_cell_mean
        y[nrow(y),ncol(y)] <- NA
        # Heat map
        chiHeat<-tibble(rep(rownames(y),times=nrow(y)),
               rep(colnames(y),each=ncol(y)),
               Mean = c(y)) %>%
          ggplot(mapping = aes(x = rep(rownames(y),times=nrow(y)),
                               y = rep(colnames(y),each=ncol(y)),
                               fill = Mean)) +
          geom_tile() +
          xlab(label = "Woman's type")+
          ylab(label= "Man's type")+
          ggtitle("Mean chi-square contribution by cell") +
          scale_fill_gradient(low = "#FFFFFF", high = "#100FCE",na.value = '#000000')
        print(chiHeat)
        
        y <- x$kl_cell_mean
        y[nrow(y),ncol(y)] <- NA
        ## Heat map
        klHeat <- tibble(rep(rownames(y),times=nrow(y)),
               rep(colnames(y),each=ncol(y)),
               Mean = c(y)) %>%
          ggplot(mapping = aes(x = rep(rownames(y),times=nrow(y)),
                               y = rep(colnames(y),each=ncol(y)),
                               fill = Mean)) +
          geom_tile() +
          xlab(label = "Woman's type")+
          ylab(label= "Man's type")+
          ggtitle("Mean KL contribution by cell") +
          scale_fill_gradient2(low = "#CE0F0F",
                                 mid = "#FFFFFF",
                                 high = "#100FCE", midpoint = 0,
                                 na.value='#000000')
	  print(klHeat)
			
}
