utils::globalVariables(c(".control.rpm"))
#' Auxiliary for Controlling rpm
#'
#' Auxiliary function as user interface for fine-tuning RPM model fitting algorithm,
#' which computes the MLPLE of the Revealed Preferences Model via optimization.
#'
#' This function is only used within a call to the \code{\link{rpm}}
#' function.
#'
#' @param init_theta vector; numeric vector of starting parameter values. This value and other possible 
#' starting values are applied to find a good optimizer. This can either have length the number of parameters
#' corresponding to the terms in the formula or in addition the equilibrium constraints.
#' @param algorithm string; The optimization algorithm to use. See \code{nloptr::nloptr.print.options()}
#' and the \code{NLopt} website for a description of the algorithms.
#' @param print_level integer; possible values: 0, 1, 2, or 3.
#' This controls how much output is shown during the
#' optimization process. Possible values: 0 (default): no output; 1:
#' show iteration number and value of objective function; 2: 1 + show
#' value of equalities/constraints; 3: 2 + show value of controls.
#' @param xtol_rel scalar; Stop when an optimization step (or an estimate of the optimum)
#' changes every parameter by less than xtol_rel multiplied by the
#' absolute value of the parameter. If there is any chance that an
#' optimal parameter is close to zero, you might want to set an absolute
#' tolerance with xtol_abs as well. Criterion is disabled if xtol_rel is
#' non-positive. Possible values: xtol_rel > 0. Default value: 1.0e-08.
#' @param ftol_rel scalar; Stop when an optimization step (or an estimate of the optimum)
#' changes the log-likelihood by less than ftol_rel multiplied by the
#' absolute value of the log-likelihood.
#' @param ftol_abs scalar; Stop when an optimization step (or an estimate of the optimum)
#' changes the log-likelihood by less than ftol_abs.
#' tolerance with xtol_abs as well. Criterion is disabled if ftol_abs is
#' non-positive. Possible values: ftol_abs > 0. Default value: 1.0e-06.
#' @param lower.bound numeric; lower bounds on the parameter estimates (that is, the beta 
#' and gamma parameters in the model). Can be a vector of the same size as the coefficient
#' vector or a single number which is used for all bounds. 
#' @param upper.bound numeric; upper bounds on the parameter estimates (that is, the beta 
#' and gamma parameters in the model). Can be a vector of the same size as the coefficient
#' vector or a single number which is used for all bounds. 
#' @param check_derivatives logical; Compare the user-supplied analytic gradients
#' with the finite difference approximations.
#' @param hessian logical; Depreciated. The negation of the `bootstrap` argument.
#' @param bootstrap logical; If  `TRUE` use a bootstrap to compute the standard errors and associated
#' covariance matrices. If `FALSE` base the standard errors and associated
#' covariance matrices on the Hessian of the
#' (constrained) log-likelihood. 
#' In all cases the extended covariance matrix is returned in \code{ext.covar.hessian}.
#' This is the matrix of parameters, log-odds of being single and the Lagrange multipliers.
#' @param maxeval integer; Stop when the number of function evaluations exceeds maxeval. This is
#' not a strict maximum: the number of function evaluations may exceed
#' maxeval slightly, depending upon the algorithm. Criterion is disabled
#' if maxeval is non-positive. Default value: 1000.
#' @param bs.maxeval integer; Stop the bootstrap optimization when the number of function evaluations exceeds bs.maxeval. This is
#' not a strict maximum: the number of function evaluations may exceed
#' bs.maxeval slightly, depending upon the algorithm. Criterion is disabled
#' if bs.maxeval is non-positive. Default value:50 
#' @param bs.xtol_rel scalar; Stop the bootstrap optimization when an optimization step (or an estimate of the optimum)
#' changes every parameter by less than bs.xtol_rel multiplied by the
#' absolute value of the parameter. See the parameter xtol_rel for details.
#' @param bs.save.data logical; Should the bootstrapped data be saved in the bootstrap return list (as components
#' \code{Xdata} and \code{Zdata}).
#' @param seed Seed value (integer) for the random number generator.  See
#' \code{\link[base]{set.seed}}
#' @param parallel.type The type of cluster to run. The typical choices are "MPI" and "PSOCK", where you
#' must have "MPI" installed to use the former. The default
#' is "PSOCK".
#' @param parallel.ncores count; Depreciated. The renamed `ncores` argument.
#' @param ncores Number of processors to use in the bootstrap computations. The default
#' is 1, that is no parallel processing.
#' @param constraints string; Additional constraints to force the proportions of singles to
#' match the (weighted) population estimates? This should not be required, but does stabilize
#' the estimates in cases where there is much uncertainty.
#' The possible values are "none" and "M_single" (the numbers of male singles of 
#' each type are reproduced). Note that adding constraints leads to 
#' over-constrained optimization which may fail.
#' @param logodds_single logical; Should the log-odds ratio of being single relative to a randomly chosen person of the same sex from the 
#' the population be returned. If FALSE the log-odds of being single relative is returned. This is a pure preference parameter.
#' @param save.data logical; Should the data be saved in the return list (as components
#' \code{Xdata} and \code{Zdata}).
#' @param robust.cov logical; Should the covariance matrix of the estimates be computed using a
#' robust method (MASS::cov.mcd)? Only use if the bootstrap is unstable.
#' @param local_opts list; list of options for \code{nloptr} sub-algorithm. See the \code{nloptr} package, but these are rarely changed.
#' @param nbootstrap integer; Number of bootstrap resamples to take in the estimation of the
#' covariance matrix of the parameter estimates.
#' @param nbootstrap.SD integer; Number of bootstrap resamples to take in the estimation of the
#' variances used in the studentized bootstrap. This is run for each nbootstrap sample and so is expensive.
#' @param large.population.bootstrap integer; If the population size exceeds \code{large.population.bootstrap} then 
#' the large population approximation is used to simulate the matchings in the bootstrap. Otherwise the 
#' small population simulation is used (including the Gale-Shapley algorithm).
#' The small population method is more accurate in smaller populations, with the default cutoff being 5000 people.
#' @param alpha proportion; Type I error rate for the confidence intervals produced by the bootstrap.
#' @return A list with arguments as components.
#' @details  Some of the arguments are not yet fully implemented.
#' It will evolve slower to incorporate more
#' arguments as the package develops.
#'
#' @seealso \code{\link{rpm}}
#' @keywords models
#' @export
control.rpm <- function(init_theta=NULL, algorithm="NLOPT_LD_SLSQP", print_level=0,
                        xtol_rel=1.0e-8, ftol_rel=1e-8, ftol_abs=1.0e-6,
                        lower.bound=-10, upper.bound=10,
                        maxeval=2000, bs.maxeval=2000, bs.xtol_rel=1.0e-8, bs.save.data=FALSE,
                        check_derivatives=FALSE, bootstrap=TRUE, hessian=FALSE, seed = NULL,
                        parallel.type="PSOCK",
                        parallel.ncores=1,
                        ncores=1,
                        constraints=c("none","M_single"),
                        logodds_single=FALSE,
                        save.data=TRUE,
                        robust.cov=FALSE,
                        local_opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-7, "maxeval"=maxeval),
                        nbootstrap=50,
                        nbootstrap.SD=20,
                        large.population.bootstrap=5000,
                        alpha=0.05
                       ) {
  formal.args <- formals(sys.function())
  if (!exists(".control.rpm")) {
    control <- list()
    for (arg in names(formal.args))
      control[arg] <- list(get(arg))
  } else{
    control <- .control.rpm
    if (!missing(algorithm)) {
      control[["init_theta"]] <- init_theta
    }
    if (!missing(algorithm)) {
      control[["algorithm"]] <- algorithm
    }
    if (!missing(print_level)) {
      control[["print_level"]] <- print_level
    }
    if (!missing(xtol_rel)) {
      control[["xtol_rel"]] <- xtol_rel
    }
    if (!missing(ftol_rel)) {
      control[["ftol_rel"]] <- ftol_rel
    }
    if (!missing(ftol_abs)) {
      control[["ftol_abs"]] <- ftol_abs
    }
    if (!missing(maxeval)) {
      control[["maxeval"]] <- maxeval
    }
    if (!missing(bs.maxeval)) {
      control[["bs.maxeval"]] <- bs.maxeval
    }
    if (!missing(bs.xtol_rel)) {
      control[["bs.xtol_rel"]] <- bs.xtol_rel
    }
    if (!missing(bs.save.data)) {
      control[["bs.save.data"]] <- bs.save.data
    }
    if (!missing(seed)) {
      control[["seed"]] <- seed
    }
    if (!missing(parallel.ncores)) {
      control[["ncores"]] <- parallel.ncores
    }
    if (!missing(ncores)) {
      control[["ncores"]] <- ncores
    }
    if (!missing(check_derivatives)) {
      control[["check_derivatives"]] <- check_derivatives
    }
    if (!missing(logodds_single)) {
      control[["logodds_single"]] <- logodds_single
    }
    if (!missing(save.data)) {
      control[["save.data"]] <- save.data
    }
    if (!missing(robust.cov)) {
      control[["robust.cov"]] <- robust.cov
    }
    if (!missing(constraints)) {
      control[["constraints"]] <- constraints
    }
    if (!missing(hessian)) {
      control[["bootstrap"]] <- hessian
    }
    if (!missing(bootstrap)) {
      control[["bootstrap"]] <- bootstrap
    }
    if (!missing(nbootstrap)) {
      control[["nbootstrap"]] <- nbootstrap
    }
    if (!missing(large.population.bootstrap )) {
      control[["large.population.bootstrap "]] <- large.population.bootstrap 
    }
    if (!missing(lower.bound)) {
      control[["lower.bound"]] <- lower.bound
    }
    if (!missing(upper.bound)) {
      control[["upper.bound"]] <- upper.bound
    }
    if (!missing(alpha)) {
      control[["alpha"]] <- alpha
    }
  }
  
  class(control) <- c("control.rpm", "control.list", "list")
  control
}
