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
#' It does this using a large-population likelihood based on ideas from Menzel (2015).
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
#' @param compare_sim string; describes which two objects are compared to compute simulated goodness-of-fit
#' statistics; valid values are \code{"sim-est"}: compares the marginal distribution of pairings in a
#' simulated sample to the \code{rpm} model estimate of the marginal distribution based on that same simulated sample;
#' \code{mod-est}: compares the marginal distribution of pairings in a
#' simulated sample to the \code{rpm} model estimate used to generate the sample
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.rpm}}, which should be consulted for specifics. 
#' @param reboot logical; if this is \code{TRUE}, the program will rerun the bootstrap at the coefficient values, rather than 
#' expect the object to contain a \code{bs.results} component with the bootstrap results run at the solution values.
#' The latter is the default for \code{rpm} fits.
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including data summary statistics.
#' @return \code{\link{gof.rpm}} returns a list consisting of the following elements:
#' \item{observed_pmf}{numeric matrix giving observed probability mass distribution over different household types}
#' \item{model_pmf}{numeric matrix giving expected probability mass distribution from \code{rpm} model}
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
#' \donttest{
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' a <- gof(fit)
#' }
#' @references Goyal, Handcock, Jackson. Rendall and Yeung (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#' Menzel, K. (2015).
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
                    empirical_p = TRUE,
                    compare_sim = 'sim-est',
                    control = object$control,
                    reboot = FALSE,
                    verbose=FALSE) 
{
  if(!is.null(control$seed)) set.seed(control$seed)
  
  model_pmf <- object$pmf_est
  observed_pmf <- object$pmf
  
  out <- list()
  out$model_pmf <- object$pmf_est
  out$observed_pmf <- object$pmf
  
  out$compare_sim <- "sim-est"
  
  # Observed Hellinger divergence
  compute_hellinger <- function(observed, expected) {
    hell <- (sqrt(observed)-sqrt(expected))^2
    return(hell)
  }
  out$obs_hellinger_cell <- compute_hellinger(observed_pmf,model_pmf)
  out$obs_hellinger <- sqrt(sum(out$obs_hellinger_cell))/sqrt(2)
  
  # Observed chi-square statistic
  compute_chi_sq <- function(N,observed,expected) {
    csc <- N*(observed-expected)^2/(expected+0.5/N)
    csc[nrow(observed),ncol(observed)] = 0
    return(csc)
  }
  out$obs_chi_sq_cell <- compute_chi_sq(object$nobs,observed_pmf,model_pmf)
  out$obs_chi_sq <- sum(out$obs_chi_sq_cell)
  
  # Observed KL divergence
  compute_kl <- function(N,observed, expected) {
    klc <- -(observed+0.5/N)*log((expected+0.5/N)/(observed+0.5/N))
    klc[nrow(observed),ncol(observed)] = 0
    klc[is.na(klc)|is.nan(klc)] <- 0
    return(klc)
  }
  out$obs_kl_cell <- compute_kl(object$nobs,observed_pmf,model_pmf)
  out$obs_kl <- sum(out$obs_kl_cell)
  
  if (empirical_p) {
    # Get it from the object bootstrap output  
    if(reboot){
    if(object$N <= control$large.population.bootstrap){
        # Use the direct (small population) method
        # These are the categories of women
        NumBetaS <- dim(object$Sd)[3]
        beta_S <- object$coefficients[1:NumBetaS]

        if(dim(object$Xd)[3]>0){
          beta_w <- object$coefficients[NumBetaS+(1:object$NumBetaW)]
        }else{
          beta_w <- NULL
        }
        if(dim(object$Zd)[3]>0){
          beta_m <- object$coefficients[NumBetaS+object$NumBetaW + (1:object$NumBetaM)]
        }else{
          beta_m <- NULL
        }

        # Use the direct (small population) method
        # These are the categories of women
        Ws <- rep(seq_along(object$pmfW),round(object$pmfW*object$num_women))
        if(length(Ws) > object$num_women+0.01){
         Ws <- Ws[-sample.int(length(Ws),size=length(Ws)-object$num_women,replace=FALSE)]
        }
        if(length(Ws) < object$num_women-0.01){
         Ws <- c(Ws,sample.int(length(object$pmfW),size=object$num_women-length(Ws),prob=object$pmfW,replace=TRUE))
        }
        Ws <- Ws[sample.int(length(Ws))]
        Ms <- rep(seq_along(object$pmfM),round(object$pmfM*object$num_men))
        if(length(Ms) > object$num_men+0.01){
         Ms <- Ms[-sample.int(length(Ms),size=length(Ms)-object$num_men,replace=FALSE)]
        }
        if(length(Ms) < object$num_men-0.01){
         Ms <- c(Ms,sample.int(length(object$pmfM),size=object$num_men-length(Ms),prob=object$pmfM,replace=TRUE))
        }
        Ms <- Ms[sample.int(length(Ms))]

        # create utility matrices
        U_star = matrix(0, nrow=length(Ws), ncol = length(Ms))
        V_star = matrix(0, nrow=length(Ms), ncol = length(Ws))
        for (ii in 1:NumBetaS) {
          U_star = U_star + object$Sd[Ws,Ms,ii] * beta_S[ii] * 0.5
        }
        if(!is.empty(beta_w)){
         for (ii in 1:object$NumBetaW) {
          U_star = U_star + object$Xd[Ws,Ms,ii] * beta_w[ii]
         }
        }
        for (ii in 1:NumBetaS) {
          V_star = V_star + t(object$Sd[Ws,Ms,ii]) * beta_S[ii] * 0.5
        }
        if(!is.empty(beta_m)){
         for (ii in 1:object$NumBetaM) {
          V_star = V_star + object$Zd[Ms,Ws,ii] * beta_m[ii]
         }
        }
  
        # adjust for outside option (remain single)

        Jw <- sqrt(object$num_men)
        Jm <- sqrt(object$num_women) 
    }

    num_Xu = nrow(object$Xu)
    num_Zu = nrow(object$Zu)
    cnW <- paste(colnames(object$Xu)[2],object$Xu[,2], sep=".")
    if(ncol(object$Xu)>2){
    for(i in 3:ncol(object$Xu)){
      cnW <- paste(cnW,paste(colnames(object$Xu)[i],object$Xu[,i], sep="."),sep='.')
    }}
    cnM <- paste(colnames(object$Zu)[2],object$Zu[,2], sep=".")
    if(ncol(object$Zu)>2){
    for(i in 2:ncol(object$Zu)){
      cnM <- paste(cnM,paste(colnames(object$Zu)[i],object$Zu[,i], sep="."),sep='.')
    }}

    if(control$ncores > 1){
      doFuture::registerDoFuture()
      future::plan(multisession, workers=control$ncores)  ## on MS Windows
      if (!is.null(object$control$seed)) {
        doRNG::registerDoRNG(object$control$seed)
      }
      if(verbose) message(sprintf("Starting parallel bootstrap using %d cores.",control$ncores))
       if(object$N > control$large.population.bootstrap){
        # Use the large population bootstrap
        bs.results <-
        foreach::foreach (b=1:control$nbootstrap, .packages=c('rpm')
        ) %dorng% {
          rpm.bootstrap.large(b, object$coefficients,
                      S=object$Sd,X=object$Xd,Z=object$Zd,sampling_design=object$sampling_design,Xdata=object$Xdata,Zdata=object$Zdata,
                      sampled=object$sampled,
                      Xid=object$Xid,Zid=object$Zid,pair_id=object$pair_id,X_w=object$X_w,Z_w=object$Z_w,
                      num_sampled=object$nobs,
                      NumBeta=object$NumBeta,NumGammaW=object$NumGammaW,NumGammaM=object$NumGammaM,
                      num_Xu=num_Xu,num_Zu=num_Zu,cnW=cnW,cnM=cnM,LB=object$LB,UB=object$UB,
                      control=control)
       }
      }else{
        # Use the small population bootstrap
        bs.results <-
        foreach::foreach (b=1:control$nbootstrap, .packages=c('rpm')
        ) %dorng% {
          rpm.bootstrap.small(b, object$coefficients,
           num_women=length(Ws), num_men=length(Ms), Jw=Jw, Jm=Jm, U_star=U_star, V_star=V_star, S=object$Sd, X=object$Xd, Z=object$Zd, pmfW=object$pmfW, pmfM=object$pmfM, object$Xu, object$Zu, num_Xu, num_Zu, cnW, cnM, object$Xid, object$Zid, object$pair_id, object$sampled, object$sampling_design, object$NumBeta, object$NumGammaW, object$NumGammaM, object$LB, object$UB, control, object$nobs)
        }
      }
      if(verbose) message(sprintf("Ended parallel bootstrap\n"))
    }else{
#     Start of the non-parallel bootstrap
      bs.results <- list()
      for(i in 1:control$nbootstrap){
        if(object$N > control$large.population.bootstrap){
        bs.results[[i]] <- rpm.bootstrap.large(b, object$coefficients,
                      S=object$Sd,X=object$Xd,Z=object$Zd,sampling_design=object$sampling_design,Xdata=object$Xdata,Zdata=object$Zdata,
                      sampled=object$sampled,
                      Xid=object$Xid,Zid=object$Zid,pair_id=object$pair_id,X_w=object$X_w,Z_w=object$Z_w,
                      num_sampled=object$nobs,
                      NumBeta=object$NumBeta,NumGammaW=object$NumGammaW,NumGammaM=object$NumGammaM,
                      num_Xu=num_Xu,num_Zu=num_Zu,cnW=cnW,cnM=cnM,LB=object$LB,UB=object$UB,
                      control=control)
        }else{
        bs.results[[i]] <- rpm.bootstrap.small(b, object$coefficients,
           num_women=length(Ws), num_men=length(Ms), Jw=Jw, Jm=Jm, U_star=U_star, V_star=V_star, S=object$Sd, X=object$Xd, Z=object$Zd, pmfW=object$pmfW, pmfM=object$pmfM, object$Xu, object$Zu, num_Xu, num_Zu, cnW, cnM, object$Xid, object$Zid, object$pair_id, object$sampled, object$sampling_design, object$NumBeta, object$NumGammaW, object$NumGammaM, object$LB, object$UB, control, object$nobs)
        }
      }
    }
    }else{
      if(is.null(object$bs.results)){
        stop("Goodness-of-Fit was called with 'reboot=FALSE', but the passed object does not contain bootstrap results. Either rerun with 'reboot=TRUE' or re-fit with 'bootstrap=TRUE'.")
      }
      bs.results <- object$bs.results
      object$bs.results <- NULL
    }
    # Get it from the object bootstrap output  
    nsim <- length(bs.results)
    sample_sizes <- unlist(lapply(bs.results, '[[', "nobs"))
    simulated_pmf_est <- array(unlist(lapply(bs.results, '[[', "pmf_est")),dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    simulated_pmf <-     array(unlist(lapply(bs.results, '[[', "pmf")),dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    
    if(verbose){
      message(sprintf("Assessing rpm model goodness-of-fit: %s",out$compare_sim))
    }
    
    # Calculate Hellinger contribution by cell for simulated pmfs
    simest_hellinger_cell <- array(NA,dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    modest_hellinger_cell <- array(NA,dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    sim_hellinger_cell <- array(NA,dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    for (it in 1:nsim) {
      simest_hellinger_cell[,,it] = compute_hellinger(simulated_pmf[,,it],simulated_pmf_est[,,it])
      modest_hellinger_cell[,,it] = compute_hellinger(simulated_pmf[,,it],model_pmf)
      sim_hellinger_cell[,,it] = simest_hellinger_cell[,,it]-modest_hellinger_cell[,,it]
    }
    
    # Calculate chi-square contribution by cell for simulated pmfs
    simest_chi_sq_cell <- array(NA, dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    modest_chi_sq_cell <- array(NA, dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    sim_chi_sq_cell <- array(NA, dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    for (it in 1:nsim) {
      simest_chi_sq_cell[,,it] = compute_chi_sq(sample_sizes[it],simulated_pmf[,,it],simulated_pmf_est[,,it])
      modest_chi_sq_cell[,,it] = compute_chi_sq(sample_sizes[it],simulated_pmf[,,it],model_pmf)
      sim_chi_sq_cell[,,it] = simest_chi_sq_cell[,,it]-modest_chi_sq_cell[,,it]
      
    }
    
    # Calculate KL contribution by cell for simulated pmfs
    simest_kl_cell <- array(NA, c(nrow(object$pmf),ncol(object$pmf),nsim))
    modest_kl_cell <- array(NA, c(nrow(object$pmf),ncol(object$pmf),nsim))
    sim_kl_cell <- array(NA, c(nrow(object$pmf),ncol(object$pmf),nsim))
    for (it in 1:nsim) {
      simest_kl_cell[,,it] <- compute_kl(sample_sizes[it],simulated_pmf[,,it],simulated_pmf_est[,,it])
      modest_kl_cell[,,it] <- compute_kl(sample_sizes[it],simulated_pmf[,,it],model_pmf)
      sim_kl_cell[,,it] <- simest_kl_cell[,,it]-modest_kl_cell[,,it]

    }
    
    # Simulated chi_sq, Hellinger and KL divergence statistics
    # This is how far the data PMF should be from the model PMF under the null
    out$hellinger_simulated <- sqrt(apply(simest_hellinger_cell, 3, sum, na.rm=TRUE))/sqrt(2)
    out$chi_sq_simulated <- apply(simest_chi_sq_cell, 3, sum, na.rm=TRUE)
    out$kl_simulated   <-   apply(simest_kl_cell,     3, sum, na.rm=TRUE)

    b <- matrix(NA,ncol=nsim,nrow=nsim)
    for(it in 1:nsim){
      for(it2 in 1:nsim){
        b[it,it2] <- sqrt(sum(compute_hellinger(simulated_pmf[,,it],simulated_pmf[,,it2]), na.rm=TRUE))/sqrt(2)
      }
    }
    diag(b) <- NA
    out$sim_hellinger_bs_pmf <- apply(b,1,mean,na.rm=TRUE)
    out$obs_hellinger_bs_pmf <- apply(simulated_pmf,3,function(x){sqrt(sum(compute_hellinger(x,object$pmf), na.rm=TRUE))/sqrt(2)})
    out$empirical_p_hellinger <- mean(apply(b,1,function(x){mean(x >= out$obs_hellinger_bs_pmf,na.rm=TRUE)}))
    h <- rep(0,nsim)
    for(i in seq_along(h)){h[i] <- mean(b >= out$obs_hellinger_bs_pmf[i],na.rm=TRUE)}
    out$empirical_p_hellinger <- mean(h)
    
    for(it in 1:nsim){
      for(it2 in 1:nsim){
        b[it,it2] <- (sum(compute_chi_sq(sample_sizes[it],simulated_pmf[,,it],simulated_pmf[,,it2]), na.rm=TRUE))
      }
    }
    diag(b) <- NA
    out$sim_chi_sq_bs_pmf <- apply(b,1,mean,na.rm=TRUE)
    out$obs_chi_sq_bs_pmf <- apply(simulated_pmf,3,function(x){(sum(compute_chi_sq(object$nobs,x,object$pmf), na.rm=TRUE))})
    out$empirical_p_chi_sq <- mean(apply(b,2,function(x){mean(x >= out$obs_chi_sq_bs_pmf,na.rm=TRUE)}))
    h <- rep(0,nsim)
    for(i in seq_along(h)){h[i] <- mean(b >= out$obs_chi_sq_bs_pmf[i],na.rm=TRUE)}
    out$empirical_p_chi_sq <- mean(h)
    
    for(it in 1:nsim){
      for(it2 in 1:nsim){
        b[it,it2] <- (sum(compute_kl(sample_sizes[it],simulated_pmf[,,it],simulated_pmf[,,it2]), na.rm=TRUE))
      }
    }
    diag(b) <- NA
    out$sim_kl_bs_pmf <- apply(b,1,mean,na.rm=TRUE)
    out$obs_kl_bs_pmf <- apply(simulated_pmf,3,function(x){(sum(compute_kl(object$nobs,x,object$pmf), na.rm=TRUE))})
    out$empirical_p_kl <- mean(apply(b,2,function(x){mean(x >= out$obs_kl_bs_pmf,na.rm=TRUE)}))
    h <- rep(0,nsim)
    for(i in seq_along(h)){h[i] <- mean(b >= out$obs_kl_bs_pmf[i],na.rm=TRUE)}
    out$empirical_p_kl <- mean(h)
    
    out$kl_simulated[is.na(out$kl_simulated)]  <- 0
    out$kl_simulated_new[is.na(out$kl_simulated_new)]  <- 0
    out$sim_chi_sq_bs_pmf[is.na(out$sim_chi_sq_bs_pmf)]  <- 0
    out$obs_chi_sq_bs_pmf[is.nan(out$obs_chi_sq_bs_pmf)]  <- 0
    out$sim_chi_sq_bs_pmf[is.nan(out$sim_chi_sq_bs_pmf)]  <- 0
    out$obs_chi_sq_bs_pmf[is.nan(out$obs_chi_sq_bs_pmf)]  <- 0
    # Remove error est from divergence statistics
    out$hellinger_simulated_new <- sqrt(apply(simest_hellinger_cell, 3, sum, na.rm=TRUE))/sqrt(2)-sqrt(apply(modest_hellinger_cell, 3, sum, na.rm=TRUE))/sqrt(2)
    out$chi_sq_simulated_new <- apply(simest_chi_sq_cell, 3, sum, na.rm=TRUE)-apply(modest_chi_sq_cell, 3, sum, na.rm=TRUE)
    out$kl_simulated_new   <-   apply(simest_kl_cell,     3, sum, na.rm=TRUE)-apply(modest_kl_cell,     3, sum, na.rm=TRUE)
    
    # Empirical p-value
    # This is rank of the observed divergence from that we should see under the null
    out$empirical_p_hellinger <- mean(out$hellinger_simulated >= out$obs_hellinger, na.rm=TRUE)
    out$empirical_p_chi_sq <- mean(out$chi_sq_simulated >= out$obs_chi_sq, na.rm=TRUE)
    out$empirical_p_kl     <- mean(out$kl_simulated >= out$obs_kl, na.rm=TRUE)
    
    # Empirical p-value
    out$empirical_p_hellinger_new <- mean(out$hellinger_simulated_new >= out$obs_hellinger, na.rm=TRUE)
    out$empirical_p_chi_sq_new <- mean(out$chi_sq_simulated_new >= out$obs_chi_sq, na.rm=TRUE)
    out$empirical_p_kl_new     <- mean(out$kl_simulated_new >= out$obs_kl, na.rm=TRUE)
    
    # Cell summaries of contributions
    # Mean, SD, Median, IQR for Hellinger
    r_dimnames <- dimnames(observed_pmf)
    out$hellinger_cell_mean <- apply(sim_hellinger_cell, 1:2, mean, na.rm=TRUE)
    dimnames(out$hellinger_cell_mean) <- r_dimnames
    
    out$hellinger_cell_sd <- apply(sim_hellinger_cell, 1:2, stats::sd, na.rm=TRUE)
    dimnames(out$hellinger_cell_sd) <- r_dimnames    
    
    out$hellinger_cell_median <- apply(sim_hellinger_cell, 1:2, stats::median, na.rm=TRUE)
    dimnames(out$hellinger_cell_median) <- r_dimnames
    
    out$hellinger_cell_iqr <- apply(sim_hellinger_cell, 1:2, stats::IQR, na.rm=TRUE)
    dimnames(out$hellinger_cell_iqr) <- r_dimnames
    
    # Mean, SD, Median, IQR for chi-sq
    out$chi_sq_cell_mean <- apply(sim_chi_sq_cell, 1:2, mean, na.rm=TRUE)
    dimnames(out$chi_sq_cell_mean) <- r_dimnames
    
    out$chi_sq_cell_sd <- apply(sim_chi_sq_cell, 1:2, stats::sd, na.rm=TRUE)
    dimnames(out$chi_sq_cell_sd) <- r_dimnames    
    
    out$chi_sq_cell_median <- apply(sim_chi_sq_cell, 1:2, stats::median, na.rm=TRUE)
    dimnames(out$chi_sq_cell_median) <- r_dimnames
    
    out$chi_sq_cell_iqr <- apply(sim_chi_sq_cell, 1:2, stats::IQR, na.rm=TRUE)
    dimnames(out$chi_sq_cell_iqr) <- r_dimnames
    
    # Mean, SD, Median, IQR for KL
    out$kl_cell_mean <- apply(sim_kl_cell, 1:2, mean, na.rm=TRUE)
    dimnames(out$kl_cell_mean) <- r_dimnames
    
    out$kl_cell_sd <- apply(sim_kl_cell, 1:2, stats::sd, na.rm=TRUE)
    dimnames(out$kl_cell_sd) <- r_dimnames
    
    out$kl_cell_median <- apply(sim_kl_cell, 1:2, stats::median, na.rm=TRUE)
    dimnames(out$kl_cell_median) <- r_dimnames
    
    out$kl_cell_iqr <- apply(sim_kl_cell, 1:2, stats::IQR, na.rm=TRUE)
    dimnames(out$kl_cell_iqr) <- r_dimnames
  }
  
  class(out) <- "gofrpm"
  return(out)
}  

#' @method print gofrpm
#' @export
print.gofrpm <- function(x, ...){
  message(sprintf("Hellinger goodness-of-fit p-value: %f", x$empirical_p_hellinger ))
  message(sprintf("Chi-Squared goodness-of-fit p-value: %f", x$empirical_p_chi_sq ))
  message(sprintf("Kullback-Liebler goodness-of-fit p-value: %f", x$empirical_p_kl ))
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
  
  hist(x$chi_sq_simulated, xlim=range(c(x$chi_sq_simulated,x$obs_chi_sq),finite=TRUE),
       probability=TRUE,
       main=expression(paste("Distribution of simulated ",chi^2)),
       xlab=expression(chi^2),nclass=20)
  lines(stats::density(x$chi_sq_simulated),col=3)
  abline(v=x$obs_chi_sq, col=2)
  
  hist(x$kl_simulated, xlim=range(c(x$kl_sim,x$obs_kl),finite=TRUE),
       probability=TRUE,
       main=expression(paste("Distribution of simulated ",KL)),
       xlab=expression(KL),nclass=20)
  lines(stats::density(x$kl_simulated),col=3)
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
  plot(chiHeat)
  
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
  plot(klHeat)
  
}
