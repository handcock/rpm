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
#' @param nsim integer; (Optional) Number of samples to draw from a large population, simulated based on the model specified
#' by \code{object}; Must be specified if \code{empirical_p} is TRUE
#' @param seed integer; (Optional) random number seed.
#' @param ncores integer; number of cores to use for parallel processing; If greater than 1
#' then the function will use parallel processing when computing empirical p value
#' @param parallel.type string; type of cluster to create. Typically, "POCK" or "MPI".
#' @param compare_sim string; describes which two objects are compared to compute simulated goodness-of-fit
#' statistics; valid values are \code{"sim-est"}: compares the marginal distribution of pairings in a
#' simulated sample to the \code{rpm} model estimate of the marginal distribution based on that same simulated sample;
#' \code{mod-est}: compares the marginal distribution of pairings in a
#' simulated sample to the \code{rpm} model estimate used to generate the sample
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
#' data(fauxmatching)
#' \donttest{
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' a <- gof(fit)
#' }
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
                    empirical_p = FALSE,
                    nsim=1000,
                    seed = NULL,
                    ncores = 1,
                    parallel.type="SOCK",
                    compare_sim = 'sim-est',
                    verbose=FALSE) 
{
  if(!is.null(seed)) set.seed(seed)
  if (empirical_p & is.null(nsim)) {
    stop("When empirical_p is TRUE, please specify the number of simulations used to calculate the empirical p value.")
  }
  
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
    csc <- N*(observed-expected)^2/(expected)
    csc[nrow(observed),ncol(observed)] = 0
    return(csc)
  }
  out$obs_chi_sq_cell <- compute_chi_sq(object$N,observed_pmf,model_pmf)
  out$obs_chi_sq <- sum(out$obs_chi_sq_cell)
  
  # Observed KL divergence
  compute_kl <- function(observed, expected) {
    klc <- -observed*log(expected/observed)
    klc[observed<1e-10] <- 0
    klc[nrow(observed),ncol(observed)] = 0
    return(klc)
  }
  out$obs_kl_cell <- compute_kl(observed_pmf,model_pmf)
  out$obs_kl <- sum(out$obs_kl_cell)
  
  if (empirical_p) {
    if (ncores > 1) {
      doFuture::registerDoFuture()
      cl <- parallel::makeCluster(ncores, type=object$control$parallel.type)
      future::plan(cluster, workers = cl)
      if(Sys.info()[["sysname"]] == "Windows"){
        future::plan(multisession)  ## on MS Windows
      }else{
       if(object$control$parallel.type!="MPI"){
        future::plan(multisession)  ## on Linux, Solaris, and macOS
       }
      }
      if (!is.null(object$control$seed)) {
        doRNG::registerDoRNG(object$control$seed)
      }
      
      pmf_sim_all <-
        foreach (b=1:nsim, .packages=c('rpm')
        ) %dorng% {
          simulated_out = microsimulate(object,
                                        pmfW_N=object$pmfW*(exp(object$gw)*object$N),
                                        pmfM_N=object$pmfM*(exp(object$gm)*object$N),
                                        large.population=FALSE)
          simulated_fit = rpm(object$formula,
                              simulated_out$population$Xdata,
                              simulated_out$population$Zdata,
                              Xid="pid",
                              Zid="pid",
                              pair_id="pair_id",
                              X_w="X_w",
                              Z_w="Z_w", 
                              pair_w="pair_w",
                              sampled="sampled",
                              sampling_design=object$sampling_design)
          return(list(countsboot=simulated_fit$counts,
                      pmfboot=simulated_fit$pmf,
                      pmfboot_est=simulated_fit$pmf_est,
                      sample_size=simulated_fit$nobs
          )
          )
        }
      parallel::stopCluster(cl)
      sample_sizes <- unlist(lapply(pmf_sim_all, '[[', 4))
      simulated_pmf_est <- array(unlist(lapply(pmf_sim_all, '[[', 3)),dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      simulated_pmf <-     array(unlist(lapply(pmf_sim_all, '[[', 2)),dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      
    } else { # not parallel
      countsboot <- array(NA,dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      sample_sizes <- rep(0, nsim)
      simulated_pmf_est <- array(NA,dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      simulated_pmf <-     array(NA,dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
      for (it in 1:nsim) {
        simulated_out = microsimulate(object,
                                      pmfW_N=object$pmfW*(exp(object$gw)*object$N),
                                      pmfM_N=object$pmfM*(exp(object$gm)*object$N),
                                      large.population=FALSE)
        simulated_fit = rpm(object$formula,
                            simulated_out$population$Xdata,
                            simulated_out$population$Zdata,
                            Xid="pid",
                            Zid="pid",
                            pair_id="pair_id",
                            X_w="X_w",
                            Z_w="Z_w", 
                            pair_w="pair_w",
                            sampled="sampled",
                            sampling_design=object$sampling_design)
        countsboot[,,it] <- simulated_fit$counts
        simulated_pmf[,,it]  <- simulated_fit$pmf
        simulated_pmf_est[,,it] <- simulated_fit$pmf_est
        sample_sizes[it] <- simulated_fit$nobs
      }
    }
    
    if(verbose){
      message(sprintf("Assessing rpm model goodness-of-fit: %s",out$compare_sim))
    }
    # Calculate Hellinger contribution by cell for simulated pmfs
    sim_hellinger_cell <- array(NA, dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    for (it in 1:nsim) {
      sim_hellinger_cell[,,it] = compute_hellinger(simulated_pmf[,,it],simulated_pmf_est[,,it])
    }
    
    # Calculate chi-square contribution by cell for simulated pmfs
    sim_chi_sq_cell <- array(NA, dim=c(nrow(object$pmf),ncol(object$pmf),nsim))
    for (it in 1:nsim) {
      sim_chi_sq_cell[,,it] = compute_chi_sq(sample_sizes[it],simulated_pmf[,,it],simulated_pmf_est[,,it])
    }
    
    # Calculate KL contribution by cell for simulated pmfs
    sim_kl_cell <- array(NA, c(nrow(object$pmf),ncol(object$pmf),nsim))
    for (it in 1:nsim) {
     #sim_kl_cell[,,it] <- compute_kl(simulated_pmf[,,it],model_pmf)
      sim_kl_cell[,,it] <- compute_kl(simulated_pmf[,,it],simulated_pmf_est[,,it])
    }
    
    # Simulated chi_sq, Hellinger and KL divergence statistics
    # This is how far the data PMF should be from the model PMF under the null
    out$hellinger_simulated <- sqrt(apply(sim_hellinger_cell, 3, sum, na.rm=TRUE))/sqrt(2)
    out$chi_sq_simulated <- apply(sim_chi_sq_cell, 3, sum, na.rm=TRUE)
    out$kl_simulated   <-   apply(sim_kl_cell,     3, sum, na.rm=TRUE)
    
    # Empirical p-value
    # This is rank of the observed divergence from that we should see under the null
    out$empirical_p_hellinger <- mean(out$hellinger_simulated >= out$obs_hellinger, na.rm=TRUE)
    out$empirical_p_chi_sq <- mean(out$chi_sq_simulated >= out$obs_chi_sq, na.rm=TRUE)
    out$empirical_p_kl     <- mean(out$kl_simulated >= out$obs_kl, na.rm=TRUE)
    
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
  print( mean(x$hellinger_sim >= x$obs_hellinger, na.rm=TRUE) )
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
