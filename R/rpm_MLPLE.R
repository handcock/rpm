#' Fit a Revealed Preference Matchings Model
#' 
#' \code{\link{rpm_MLPLE}} estimates the parameters of a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using a large-population approximation to the likelihood based on ideas from Menzel (2015).
#' 
#' It is usually called via the \code{\link{rpm}} function.
#' @param formula formula; an \code{\link{formula}} object, of the form \code{
#' ~ <model terms>}. For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param Xdata data.frame for women. Each row is a woman, each column is a variable on that women
#' or her partnerships. It must contain the women's ID variable (see \code{Xid}) and
#' a variable with the ID of the women's partner. If the women is single the men's ID should be NA.
#' @param Zdata data.frame for men. Each row is a man, each column is a variable on that men
#' It must contain the men's ID variable (see \code{Zid}).
#' @param Xid string The name of the variable in \code{Xdata} containing the IDs of the women.
#' @param Zid string The name of the variable in \code{Zdata} containing the IDs of the men.
#' @param pair_id string The name of the variable in \code{Xdata} containing the ID of the men paired with the women in
#' \code{Xid}. If the women is not paired it must be NA.
#' @param X_w string The name of the variable in \code{Xdata} containing the individual weight of the women.
#' If this is NULL then it is assumed the sample is unweighted from a population with 2000 women in it.
#' @param Z_w string The name of the variable in \code{Zdata} containing the individual weight of the man
#' If this is NULL then it is assumed the sample is unweighted from a population with 2000 men in it.
#' @param pair_w string The name of the variable in \code{Xdata} containing the pair weight of that women.
#' If the women is not paired it should be NA.
#' If this is NULL then it is computed from the individual weights using the \code{sampling_design}.
#' Note that the pair weights currently do not play a role in the estimation. They do in the quasi-likelihood version of the code.
##' If this is NULL then it is assumed the sample is unweighted from a population with 2000 men in it.
#' @param sampled string The name of the logical variable in \code{Xdata} and \code{Zdata} containing the 
#' indicator that the person was sampled directly (as distinct from being included as the match of a directly sampled
#' person. All single people are directly sampled. 
#' @param sampling_design string; The name of the sampling protocol used to select the survey data. Valid values are
#' \code{"stock-flow"} (default) (individuals are sampled, data contains both
#' singles and couples);
#' \code{"stock-stock"} (households are sampled, each household can be a single or a couple);
#' \code{"census"} (the sample is a census of the population of people).
## (only couples are included in the data).
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.rpm}}, which should be consulted for specifics. 
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including data summary statistics.
#' @details The pairings are determined by the \code{pair_id} variable in \code{Xdata}. 
#' If that variable is NA then the women is
#' assumed to be single. If men are listed in \code{Zdata} and are not partnered then they are assumed single.
#' Weights are specified by three optional variables in \code{Xdata}.
#' \itemize{
#' \item{X_w}{This is character string of the name of the weight variable for women. The sum of the weights should be the
#' number of women in the population.}
#' \item{Z_w}{This is character string of the name of the weight variable for men. The sum of the weights should be the
#' number of men in the population.}
#' \item{pair_w}{This is character string of the name of the weight variable for pairs.}
#' }
#' @return \code{\link{rpm}} returns an object of class \code{\link{rpm.object}}
#' that is a list consisting of the following elements: 
#' \item{coefficients}{The bias-corrected bootstrap estimate of \eqn{\theta}, the vector of
#' coefficients for the model parameters. This includes the model \eqn{\beta} and the model \eqn{\Gamma}.}
#' \item{loglik}{The value of the maximized log-likelihood.}
#' \item{exitflag}{integer value with the status of the optimization (4 is success as 
#' \code{xtol_rel} or \code{xtol_abs} was reached). Other codes are 1 = generic success; 2 = optimization stopped because 
#' \code{ftol_rel} or \code{ftol_abs} was reached; 3 = optimization stopped 
#' because \code{stopval} was reached; 4 = optimization stopped because \code{xtol_rel} or \code{xtol_abs} was reached; 
#' 5 = optimization stopped because
#' \code{maxeval} was reached; 6 = optimization stopped because \code{maxtime} was reached.}
#' \item{call}{the call that was made to \code{nloptr}.}
#' \item{x0}{vector with starting values for the optimization.}
#' \item{message}{more informative message with the status of the optimization.}
#' \item{iterations}{number of iterations that were executed.}
#' \item{objective}{value if the objective function in the solution.}
#' \item{solution}{optimal value of the controls.}
#' \item{version}{version of NLopt that was used.}
#' \item{covar}{Approximate covariance matrix of the estimates.}
#' \item{eq}{Values from the equality constraints. Larger values indicate non-convergence.}
#' @seealso control.rpm, summary.rpm, print.rpm
#' @references Goyal, Handcock, Jackson. Rendall and Yeung (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#' Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' library(rpm)
#' \donttest{
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled",sampling_design="stock-flow")
#' summary(fit)
#' }
#' @importFrom MASS cov.mcd
#' @export rpm_MLPLE
rpm_MLPLE <- function(formula, Xdata, Zdata,
    Xid=NULL, Zid=NULL, pair_id=NULL,
    X_w=NULL, Z_w=NULL, pair_w=NULL,
    sampled=NULL,
    sampling_design="stock-flow",
    control=control.rpm(), verbose=FALSE){

    if(!is.null(control$seed)) set.seed(control$seed)

    if(is.null(Xid)){
      warning("The variable name in Xdata of the women's IDs is missing and is presumed to be 'Xid'.")
      Xid <- "Xid"
      if(is.null(Xdata[["Xid"]])){
        stop("The variable name in Xdata of the women's IDs must be specified as the character string 'Xid'.")
      }
    }
    if(is.null(Zid)){
      warning("The variable name in Zdata of the men's IDs is missing and is presumed to be 'Zid'.")
      Zid <- "Zid"
      if(is.null(Zdata[["Zid"]])){
        stop("The variable name in Zdata of the women's IDs must be specified as the character string 'Xid'.")
      }
    }
    if(is.null(pair_id)){
      warning("The variable name in Xdata of the men's IDs of the women's partners has not been specified. It is presumed to be 'pair_id'.")
      pair_id <- "pair_id"
      if(is.null(Xdata[["pair_id"]])){
        stop("The variable 'pair_id' in Xdata of the men's IDs of the women's partners is NULL.")
      }
    }
    if(is.null(sampled) & sampling_design != "census"){
      stop("The variable name in Xdata and Zdata of the (directly) sampled indicator must be specified as the character string 'sampled'.")
    }
    if(is.null(X_w) & sampling_design != "census"){
      warning("The variable name in Xdata of the women's weights has not been specified. It is presumed to be 'X_w'.")
      X_w <- "X_w"
      if(is.null(Xdata[["X_w"]])){
    #   This is the unweighted case
        warning("The variable name in Xdata of the women's weights has not been specified. It is specified by the character string 'X_w'. We will presume the sample is unweighted from a population with 2000 women in it.")
        nw <- 2000
        Xdata[["X_w"]]=rep(nw/nrow(Xdata),length=nrow(Xdata))
      }
    }

    if(is.null(Z_w) & sampling_design != "census"){
      warning("The variable name in Zdata of the men's weights has not been specified. It is presumed to be 'Z_w'.")
      Z_w <- "Z_w"
      if(is.null(Zdata[["Z_w"]])){
      # This is the unweighted case
        warning("The variable name in Zdata of the men's weights has not been specified. It is specified by the character string 'Z_w'. We will presume the sample is unweighted from a population with 2000 men in it.")
        nm <- 2000
        Zdata[["Z_w"]]=rep(nm/nrow(Zdata),length=nrow(Zdata))
      }
    }

    if(is.null(X_w)){ X_w <- "X_w" }
    if(is.null(Z_w)){ Z_w <- "Z_w" }
    
    if(sampling_design == "stock-flow" && all(Xdata[,sampled]) && all(Zdata[,sampled])){
      warning("The sampling design is specified as stock-flow, but all men and women are designated as directly sampled.")
    }

    # coerce IDs to character to aid comparison
    Xdata[,Xid] <- as.character(Xdata[,Xid])
    Zdata[,Zid] <- as.character(Zdata[,Zid])
    Xdata[,pair_id] <- as.character(Xdata[,pair_id])
    Zdata[,pair_id] <- as.character(Zdata[,pair_id])

    M_paired_to_W <- Zdata[!is.na(Zdata[,pair_id]),pair_id] %in%  Xdata[,Xid]
    W_paired_to_M <- Xdata[!is.na(Xdata[,pair_id]),pair_id] %in%  Zdata[,Zid]
    W_pair_id <- Xdata[!is.na(Xdata[,pair_id]),pair_id]
    W_duplicated_pair_id <- W_pair_id[duplicated(W_pair_id)]
    M_pair_id <- Zdata[!is.na(Zdata[,pair_id]),pair_id]
    M_duplicated_pair_id <- M_pair_id[duplicated(M_pair_id)]

    if(!all(W_paired_to_M)){
      message("The women's data, Xdata, contains people with pair IDs that do not appear in the men's data, Zdata. These are:")
      a <- Xdata[!is.na(Xdata[,pair_id]),][!W_paired_to_M,]
      message_print(a)
      message("Please correct this before rerunning.")
    }

    if(!all(M_paired_to_W)){
      message("The men's data, Zdata, contains people with pair IDs that do not appear in the women's data, Xdata. These are:")
      a <- Zdata[!is.na(Zdata[,pair_id]),][!M_paired_to_W,]
      message_print(a)
      message("Please correct this before rerunning.")
    }

    if(!is.empty(W_duplicated_pair_id)){
      message("The women's data, Xdata, contains duplicated pair IDs. These are:")
      a <- Xdata[!is.na(Xdata[,pair_id]),][W_pair_id %in% W_duplicated_pair_id,]
      a <- a[order(a[,pair_id]),]
      message_print(a)
      message("Please correct this before reruning.")
    }
    if(!is.empty(M_duplicated_pair_id)){
      message("The men's data, Zdata, contains duplicated pair IDs. These are:")
      a <- Zdata[!is.na(Zdata[,pair_id]),][M_pair_id %in% M_duplicated_pair_id,]
      a <- a[order(a[,pair_id]),]
      message_print(a)
      message("Please correct this before reruning.")
    }
    if(!is.empty(W_duplicated_pair_id) | !is.empty(M_duplicated_pair_id) | !all(W_paired_to_M) | !all(M_paired_to_W)){
      stop()
    }

    # 1) parse the formula
    # intercept is always added as the first column
    if(is.na(match("Int", colnames(Xdata)[1]))){
      Xdata <- cbind(1, Xdata)
      colnames(Xdata)[1] <- "Int"
    }
    if(is.na(match("Int", colnames(Zdata)[1]))){
      Zdata <- cbind(1, Zdata)
      colnames(Zdata)[1] <- "Int"
    }
    
    model.terms<-list_rhs.formula(formula)
    
    # 2) get the subset of the variables that are relevant according to the formula
    # Define the variables used in the model (and hence the unique classes of partners)
    # This is typically a subset of the available variables
    model_vars_fn <- function(x){
     if(length(x)>1){
       as.character(x[[2]])
     }else{
       NULL
     }
    }
    model_vars <- c("Int", unlist(unique(lapply(model.terms, model_vars_fn))))
    
    W_matched_vars <- model_vars %in% colnames(Xdata)    
    if(!all(W_matched_vars)){
      message("Some of the variables in the formula do not appear in the women's data, Xdata. These are:")
      message(model_vars[!W_matched_vars])
      message("Please correct this before reruning.")
    }
    M_matched_vars <- model_vars %in% colnames(Zdata)    
    if(!all(M_matched_vars)){
      message("Some of the variables in the formula do not appear in the men's data, Zdata. These are:")
      message(model_vars[!M_matched_vars])
      message("Please correct this before reruning.")
    }
    if(!all(W_matched_vars) | !all(W_matched_vars)){
      stop()
    }

    # 3) Compute the marginal distributions of women's and men's types
    Xu <- unique(Xdata[,model_vars])
    if(length(Xu)==1){
      Xu <- as.matrix(Xu[do.call(order, as.data.frame(Xu)),drop=FALSE])
    }else{
      Xu <- Xu[do.call(order, as.data.frame(Xu)),]
    }
    Zu <- unique(Zdata[,model_vars])
    if(length(Zu)==1){
      Zu <- as.matrix(Zu[do.call(order, as.data.frame(Zu)),drop=FALSE])
    }else{
      Zu <- Zu[do.call(order, as.data.frame(Zu)),]
    }

    if(ncol(Xu)>1){
      cnW <- paste(colnames(Xu)[2],Xu[,2], sep=".")
      if(ncol(Xu)>2){
      for(i in 3:ncol(Xu)){
        cnW <- paste(cnW,paste(colnames(Xu)[i],Xu[,i], sep="."),sep='.')
      }}
    }else{cnW <- "Int"}
    if(ncol(Zu)>1){
      cnM <- paste(colnames(Zu)[2],Zu[,2], sep=".")
      if(ncol(Zu)>2){
      for(i in 2:ncol(Zu)){
        cnM <- paste(cnM,paste(colnames(Zu)[i],Zu[,i], sep="."),sep='.')
    }}
    }else{cnM <- "Int"}
    
    # 4) Create joint PMF
    # Xtype: group membership for women (one for each woman in the pop)
    d <- matrix(match(t(Xdata[,model_vars]),t(Xu)),ncol=ncol(Xu),byrow=TRUE)
    d <- t(apply(d,1,function(x){x/seq_along(model_vars)}))
    Xtype <- matrix(match(t(Xu),t(Xu)),ncol=ncol(Xu),byrow=TRUE)
    Xtype <- t(apply(Xtype,1,function(x){x/seq_along(model_vars)}))
    if(length(model_vars)>1){
      Xtype <- match(data.frame(t(d)),data.frame(t(Xtype)))
    }else{
      Xtype <- as.vector(d)
    }
    Xdata[["Xtype"]] <- Xtype
    rm(Xtype)
    # Ztype: group membership for men (one for each man in the pop)
    d <- matrix(match(t(Zdata[,model_vars]),t(Zu)),ncol=ncol(Zu),byrow=TRUE)
    d <- t(apply(d,1,function(x){x/seq_along(model_vars)}))
    Ztype <- matrix(match(t(Zu),t(Zu)),ncol=ncol(Zu),byrow=TRUE)
    Ztype <- t(apply(Ztype,1,function(x){x/seq_along(model_vars)}))
    if(length(model_vars)>1){
      Ztype <- match(data.frame(t(d)),data.frame(t(Ztype)))
    }else{
      Ztype <- as.vector(d)
    }
    Zdata[["Ztype"]] <- Ztype
    rm(Ztype)
    
    if (!is.null(X_w) & sampling_design != "census") {
        Xdata$X_w = Xdata[,X_w]
        Xdata$sampled = Xdata[,sampled]
    }
    if (!is.null(Z_w) & sampling_design != "census") {
        Zdata$Z_w = Zdata[,Z_w]
        Zdata$sampled = Zdata[,sampled]
    }

    num_Xu = nrow(Xu)
    num_Zu = nrow(Zu)

    mc <- rpm_make_counts(Xdata, Zdata, sampling_design, sampled, Xid, Zid, pair_id, X_w, Z_w, Xu, Zu, verbose)
    pmf <- mc$pmf
    num_sampled=mc$num_sampled
    counts=mc$counts
    pmfW=mc$pmfW
    pmfM=mc$pmfM
    pmfN=mc$pmfN
    N=mc$N
    gw=mc$gw
    gm=mc$gm
    num_women=mc$num_women
    num_men=mc$num_men

    # 5) create model matrix
    modelmat <- rpm.model.matrix(model.terms, Xu, Zu)
    
    S <- modelmat$S
    X <- modelmat$X
    Z <- modelmat$Z

    NumBetaS <- dim(S)[3]
    NumBetaW <- dim(X)[3]
    NumBetaM <- dim(Z)[3]
    NumBeta <- NumBetaS + NumBetaW + NumBetaM
    NumGammaW <- num_Xu
    NumGammaM <- num_Zu
    NumGamma <- NumGammaW + NumGammaM
    
    # 6) Set init_theta to a good starting value to save time
    init_theta = control[["init_theta"]]
    if(is.null(init_theta)|!is.vector(init_theta, mode = "numeric")|!(length(init_theta) %in% c(NumBeta,NumBeta+NumGamma))){
        init_theta <- 
                       c(log(sum(pmf[-nrow(pmf),-ncol(pmf)])/(1-sum(pmf[-nrow(pmf),-ncol(pmf)]))),
                        rep(0,NumBeta-1),
                       -log(c(2*apply(pmf[-nrow(pmf),-ncol(pmf),drop=FALSE],1,sum)/pmf[-nrow(pmf),ncol(pmf),drop=FALSE],
                              2*apply(pmf[-nrow(pmf),-ncol(pmf),drop=FALSE],2,sum)/pmf[nrow(pmf),-ncol(pmf),drop=FALSE])))
    }else{
     tmp <- init_theta
     if(identical(length(init_theta),NumBeta+NumGamma)){
       init_theta <- tmp
     }else{if(identical(length(init_theta),NumBeta)){
       init_theta <- c(tmp,
                       -log(c(2*apply(pmf[-nrow(pmf),-ncol(pmf),drop=FALSE],1,sum)/pmf[-nrow(pmf),ncol(pmf),drop=FALSE],
                              2*apply(pmf[-nrow(pmf),-ncol(pmf)],2,sum)/pmf[nrow(pmf),-ncol(pmf)])))
     }}
    }
    init_theta[is.na(init_theta) | is.nan(init_theta)] <- 0
    nstr = c(modelmat$Snames, modelmat$Xnames, modelmat$Znames)
    names(init_theta) <- c(nstr, paste0("LOD_Single.W.",cnW), paste0("LOD_Single.M.",cnM))
    if(verbose){
      message("Initial estimate set to:")
      message_print(init_theta)
    }

    # core algorithm begins
    
    constraints <- control[["constraints"]]
    constraints <- match.arg(control[["constraints"]], c("none","M_single"))
    constraints <- match(constraints, c("M_single","none")) - 1

    if(verbose){
      message(switch(constraints+1,
        sprintf("Constraints: Adding single men's proportions by type of man."),
        sprintf("No additional constraints applied.")
       )
      )
    }

    if(length(control$lower.bound)<(NumBeta+NumGamma)){
      LB <- c(rep(control$lower.bound[1],NumBeta),rep(control$lower.bound[1],NumGamma))
    }else{
      LB <- control$lower.bound[1:(NumBeta+NumGamma)]
    }
    if(length(control$upper.bound)<(NumBeta+NumGamma)){
      UB <- c(rep(control$upper.bound[1],NumBeta),rep(control$upper.bound[1],NumGamma))
    }else{
      UB <- control$upper.bound[1:(NumBeta+NumGamma)]
    }
    if(any(LB>UB) | all(LB==UB)){
      stop("The lower bounds on the estimates should be at most the upper bounds. Please respecify them.")
    }

    loglikfun=loglikfun_nog
    gloglikfun=gloglikfun_nog
    eqfun=eqfun_nog
    jeqfun=jeqfun_nog

    test_eqcond=eqfun(init_theta, Sd=S, Xd=X, Zd=Z, NumGammaW, NumGammaM, pmfW, pmfM, pmf, pmfN, gw, gm, N, sampling_design, constraints)

    init_theta[init_theta>=UB] = UB[init_theta>=UB]-0.01
    init_theta[init_theta<=LB] = LB[init_theta<=LB]+0.01
    out.text <- capture.output(
        out.list <- nloptr::nloptr(x0=init_theta, eval_f=loglikfun, 
                  eval_grad_f=gloglikfun,
                  eval_g_eq=eqfun,  eval_jac_g_eq=jeqfun,
                  lb=LB,ub=UB,
                  Sd=S,Xd=X,Zd=Z,NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                  pmfW=pmfW, pmfM=pmfM, pmf=pmf, counts=pmfN, gw=gw, gm=gm, N=N,
                  sampling=sampling_design, constraints=constraints,
                  opts=control)
    )
    if(any(startsWith("Error",c(" ",out.text)))){
      message(sprintf("Optimization is overly constrained. Estimates may be unstable."))
    }

    names(out.list$solution) <- names(init_theta)
    
    if(verbose){
      message("Maximized log-likelihood:")
      message(sprintf("log-likelihood = %s",
               format(-out.list$objective,nsmall=5,digits=6)))
      a <- out.list$solution
      names(a) <- names(init_theta)
      message_print(a)
    }

    out.hessian <- rpm.hessian_nog(out.list$solution,Sd=S,Xd=X,Zd=Z,
             NumBeta=NumBeta, NumGamma=NumGamma,
             NumGammaW=NumGammaW, NumGammaM=NumGammaM,
             pmfW=pmfW, pmfM=pmfM, pmf=pmf, counts=pmfN, gw=gw, gm=gm, N=N,
             sampling=sampling_design, constraints=constraints,verbose=verbose)

    bs.results <- NULL
    if(!control$bootstrap){
     covar.val <- function(out.hessian){
      ifelse(all(!is.na(diag(out.hessian)))&all(diag(out.hessian)>0),0.01*sum(log(diag(out.hessian))),1000000000)
     }
     llik.hessian <- function(out,out.hessian){
       ifelse(!is.null(out) & out$status >= 1 & out$status <=4 & out$objective > 0,out$objective,1000000000)
     }
     lliks <- llik.hessian(out.list,out.hessian$covar)
     covars <- covar.val(out.hessian$covar)
     out <- out.list
     out$covar <- out.hessian$covar
     out$ext.covar <- out.hessian$ext.covar
     out$covar.unconstrained <- out.hessian$covar.unconstrained
     out$ext.covar.hessian <- out.hessian$ext.covar
     
    }else{
# no hessian version
     llik.nohessian <- function(out){
       ifelse(!is.null(out) && (out$status >= 1 & out$status <=4 & out$objective > 0),out$objective,1000000000)
     }
     lliks <- llik.nohessian(out.list)
     out <- out.list
     out$covar <- out.hessian$covar
     out$ext.covar <- out.hessian$ext.covar
     out$covar.unconstrained <- out.hessian$covar.unconstrained
     out$ext.covar.hessian <- out.hessian$ext.covar
    }
# end hessian covar

    th_hat <- out.list$solution

    if(inherits(out,"try-error") || any(is.na(out$solution))){
      message("MLPLE failed; restarting at different starting values.\n")
      if(inherits(out,"try-error")){
         message("MLPLE still contains NAs.\n")
      }
    }
  
    out$loglik <- -out$objective
    out$exitflag <- out$status
    if(verbose){ 
      message("eq values:")
      a <- eqfun(th_hat, Sd=S, Xd=X, Zd=Z, NumGammaW, NumGammaM, pmfW, pmfM, pmf, pmfN, gw, gm, N, sampling_design,constraints)
      message(a)
      message(round(th_hat,2))
    }
        
    if(is.nan(out$loglik) | is.infinite(out$loglik)){
     warning("The fit did not find a valid starting value and the estimates are not meaningful. The likelihood is likely over-constrained. If possible, try setting control.rpm(constraints='none').")
    }

    out$eq = eqfun(th_hat, Sd=S, Xd=X, Zd=Z, NumGammaW, NumGammaM, pmfW, pmfM, pmf, pmfN, gw, gm, N, sampling_design,constraints)

    pmf_est <- exp(augpmfnew(beta=th_hat[1:NumBeta],
                GammaW=th_hat[NumBeta+(1:NumGammaW)], 
                GammaM=th_hat[(NumBeta+NumGammaW)+(1:NumGammaM)],
                S=S, X=X, Z=Z,
                pmfW=pmfW, pmfM=pmfM, gw=gw, gm=gm))
    pmf_est[nrow(pmf_est),ncol(pmf_est)] <- 0

    PMF_SW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
    out$PMF_SW <- PMF_SW / (PMF_SW + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],1,sum))
    PMF_SM <- pmf_est[ nrow(pmf_est),-ncol(pmf_est),drop=FALSE]
    out$PMF_SM <- PMF_SM / (PMF_SM + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],2,sum))
    PMF_PW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
    out$PMF_PW <- 1 - out$PMF_SW
    out$PMF_PM <- 1 - out$PMF_SM 

    pmf_est[-nrow(pmf_est), -ncol(pmf_est)] <- 2*pmf_est[-nrow(pmf_est), -ncol(pmf_est)]
    pmf_est <- pmf_est/sum(pmf_est)

    names(out$PMF_SW) <- paste0("PMF_Single.W.",cnW)
    names(out$PMF_SM) <- paste0("PMF_Single.M.",cnM)
    names(out$PMF_PW) <- paste0("PMF_Partnered.W.",cnW)
    names(out$PMF_PM) <- paste0("PMF_Partnered.M.",cnM)
    if(control$logodds_single){
      b <- sum(out$PMF_SW*pmfW)
      out$LOGODDS_SW <- th_hat[(NumBeta+1):(NumBeta+NumGammaW)] - log(b/(1-b))
      b <- sum(out$PMF_SM*pmfM)
      out$LOGODDS_SM <- th_hat[(NumBeta+NumGammaW+1):(NumBeta+NumGammaW+NumGammaM)] - log(b/(1-b))
    }else{
      out$LOGODDS_SW <- log(out$PMF_SW/(1-out$PMF_SW))
      out$LOGODDS_SM <- log(out$PMF_SM/(1-out$PMF_SM))
    }
    names(out$LOGODDS_SW) <- paste0("LOD_Single.W.",cnW)
    names(out$LOGODDS_SM) <- paste0("LOD_Single.M.",cnM)

    out$solution[NumBeta+(1:NumGammaW)] <- out$LOGODDS_SW
    out$solution[NumBeta+NumGammaW+(1:NumGammaM)] <- out$LOGODDS_SM

    pmfN_households <- pmfN
    pmfN_households[-nrow(pmfN), -ncol(pmfN)] <- 0.5*pmfN[-nrow(pmfN), -ncol(pmfN)]
    pmf_households <- pmfN_households / sum(pmfN_households)
    out$loglik <- stats::dmultinom(x=pmfN_households,prob=pmf_est,log=TRUE)

    if(formula != stats::as.formula(~1)){
      PMF_SW <- pmf_households[-nrow(pmfN), ncol(pmfN)]
      PMF_SM <- pmf_households[ nrow(pmfN),-ncol(pmfN)]
      a <- outer(PMF_SW,PMF_SM,"*")
      a <- a / sum(a)
      pmf_null <- pmf_households
      pmf_null[-nrow(pmfN),-ncol(pmfN)] <- a*sum(pmf_null[-nrow(pmfN),-ncol(pmfN)])
      out$loglik.null <- stats::dmultinom(x=pmfN_households,prob=pmf_null,log=TRUE)
    }else{
      out$loglik.null <- out$loglik
      out$null_coefficients = out$coefficients     
    }

    out$counts=counts
    out$Sd=S; out$Xd=X; out$Zd=Z
    out$control <- control
    out$nobs <- num_sampled

    out$formula <- formula
    out$sampled <- sampled
    out$Xid <- Xid
    out$Zid <- Zid
    out$pair_id <- pair_id
    out$X_w <- X_w
    out$Z_w <- Z_w
    out$Xu <- Xu
    out$Zu <- Zu
    out$gw <- gw
    out$gm <- gm
    out$N <- N
    out$LB <- LB
    out$UB <- UB
    out$sampling_design <- sampling_design
    
    if (!control$bootstrap) {
     out$analyticalCI <- cbind(out$solution[1:NumBeta]-1.96*sqrt(diag(out$covar)[1:NumBeta]),
                               out$solution[1:NumBeta]+1.96*sqrt(diag(out$covar)[1:NumBeta]))
     colnames(out$analyticalCI) = c('LB','UB')
    }

    if(control$bootstrap){
     # so use the bootstrap
     # First save the hessian-based one
     out.hessian <- rpm.hessian_nog(th_hat,Sd=S,Xd=X,Zd=Z,
                  NumBeta=NumBeta, NumGamma=NumGamma,
                  NumGammaW=NumGammaW, NumGammaM=NumGammaM,
                  pmfW=pmfW, pmfM=pmfM, pmf=pmf, counts=pmfN, gw=gw, gm=gm, N=N,
                  sampling=sampling_design, constraints=constraints,verbose=verbose)
     out$covar <- out.hessian$covar
     out$ext.covar <- out.hessian$ext.covar
     out$covar.unconstrained <- out.hessian$covar.unconstrained
     out$ext.covar.hessian <- out.hessian$ext.covar

     if(control$ncores > 1){
      doFuture::registerDoFuture()
      if(Sys.info()[["sysname"]] == "Windows"){
        future::plan(multisession, workers=control$ncores)  ## on MS Windows
      }else{
        if(control$parallel.type!="MPI"){
          future::plan(multisession, workers=control$ncores)     ## on Linux, Solaris, and macOS
        }
      }
      if (!is.null(control$seed)) {
        doRNG::registerDoRNG(control$seed)
      }
    }

    if (control$nbootstrap<length(c(out$solution[1:NumBeta],out$LOGODDS_SW,out$LOGODDS_SM))){
      if (control$nbootstrap>0){
        control$nbootstrap = round(length(c(out$solution[1:NumBeta],out$LOGODDS_SW,out$LOGODDS_SM))+1/2*length(c(out$solution[1:NumBeta],out$LOGODDS_SW,out$LOGODDS_SM)))
        message(sprintf("nbootstrap changed to %d so that it is not too small to be meaningful.\n",control$nbootstrap))
      }else{
        control$nbootstrap = 2
      }
    }
    out.boot <- matrix(0,ncol=length(out$solution),nrow=control$nbootstrap)
    out.boot.LOGODDS <- matrix(0,ncol=length(c(out$LOGODDS_SW,out$LOGODDS_SM)),nrow=control$nbootstrap)
    out.boot_SD <- matrix(0,ncol=length(c(out$solution[1:NumBeta],out$LOGODDS_SW,out$LOGODDS_SM)),nrow=control$nbootstrap)
	tstar.boot <- matrix(0,ncol=NumBeta,nrow=control$nbootstrap)
    colnames(out.boot) <- names(out$solution)
    colnames(tstar.boot) <- names(out$solution)[1:NumBeta]

    if(N <= control$large.population.bootstrap){
        # Use the direct (small population) method
        # These are the categories of women
        beta_S <- out$solution[1:NumBetaS]

        if(dim(X)[3]>0){
          beta_w <- out$solution[NumBetaS+(1:NumBetaW)]
        }else{
          beta_w <- NULL
        }
        if(dim(Z)[3]>0){
          beta_m <- out$solution[NumBetaS+NumBetaW + (1:NumBetaM)]
        }else{
          beta_m <- NULL
        }

        # Use the direct (small population) method
        # These are the categories of women
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
        for (ii in 1:NumBetaS) {
          U_star = U_star + S[Ws,Ms,ii] * beta_S[ii] * 0.5
        }
        if(!is.empty(beta_w)){
         for (ii in 1:NumBetaW) {
          U_star = U_star + X[Ws,Ms,ii] * beta_w[ii]
         }
        }
        for (ii in 1:NumBetaS) {
          V_star = V_star + t(S[Ws,Ms,ii]) * beta_S[ii] * 0.5
        }
        if(!is.empty(beta_m)){
         for (ii in 1:NumBetaM) {
          V_star = V_star + Z[Ms,Ws,ii] * beta_m[ii]
         }
        }
  
        # adjust for outside option (remain single)
        Jw <- sqrt(num_men)
        Jm <- sqrt(num_women) 
    }

    if(control$ncores > 1){
      if(verbose) message(sprintf("Starting parallel bootstrap using %d cores.",control$ncores))
       if(N > control$large.population.bootstrap){
        # Use the large population bootstrap
        bs.results <-
        foreach::foreach (b=1:control$nbootstrap, .packages=c('rpm')
        ) %dorng% {
        rpm.bootstrap.large(b, out$solution,
                      S=S,X=X,Z=Z,sampling_design=sampling_design,Xdata=Xdata,Zdata=Zdata,
                      sampled=sampled,
                      Xid=Xid,Zid=Zid,pair_id=pair_id,X_w=X_w,Z_w=Z_w,
                      num_sampled=num_sampled,
                      NumBeta=NumBeta,NumGammaW=NumGammaW,NumGammaM=NumGammaM,
                      num_Xu=num_Xu,num_Zu=num_Zu,cnW=cnW,cnM=cnM,LB=LB,UB=UB,
                      control=control)
       }
      }else{
        bs.results <-
        foreach::foreach (b=1:control$nbootstrap, .packages=c('rpm')
        ) %dorng% {
        rpm.bootstrap.small(b, out$solution,
           num_women=length(Ws), num_men=length(Ms), Jw=Jw, Jm=Jm, U_star=U_star, V_star=V_star, S=S, X=X, Z=Z, pmfW=pmfW, pmfM=pmfM, Xu, Zu, num_Xu, num_Zu, cnW, cnM, Xid, Zid, pair_id, sampled, sampling_design, NumBeta, NumGammaW, NumGammaM, LB, UB, control, num_sampled)
        }
      }
      for(i in 1:control$nbootstrap){
        out.boot[i,] <- bs.results[[i]]$est
        out.boot.LOGODDS[i,1:NumGammaW] <- bs.results[[i]]$LOGODDS_SW
        out.boot.LOGODDS[i,NumGammaW+(1:NumGammaM)] <- bs.results[[i]]$LOGODDS_SM
      }
      if(verbose) message(sprintf("Ended parallel bootstrap\n"))
    }else{
#     Start of the non-parallel bootstrap
      bs.results <- list()
      for(i in 1:control$nbootstrap){
        if(N > control$large.population.bootstrap){
         bs.results[[i]] <- rpm.bootstrap.large(1, out$solution,
                      S=S,X=X,Z=Z,sampling_design=sampling_design,Xdata=Xdata,Zdata=Zdata,
                      sampled=sampled,
                      Xid=Xid,Zid=Zid,pair_id=pair_id,X_w=X_w,Z_w=Z_w,
                      num_sampled=num_sampled,
                      NumBeta=NumBeta,NumGammaW=NumGammaW,NumGammaM=NumGammaM,
                      num_Xu=num_Xu,num_Zu=num_Zu,cnW=cnW,cnM=cnM,LB=LB,UB=UB,
                      control=control)
        }else{
         bs.results[[i]] <- rpm.bootstrap.small(1, out$solution,
           num_women=length(Ws), num_men=length(Ms), Jw=Jw, Jm=Jm, U_star=U_star, V_star=V_star, S=S, X=X, Z=Z, pmfW=pmfW, pmfM=pmfM, Xu, Zu, num_Xu, num_Zu, cnW, cnM, Xid, Zid, pair_id, sampled, sampling_design, NumBeta, NumGammaW,NumGammaM, LB, UB, control, num_sampled)
        }
        out.boot[i,] <- bs.results[[i]]$est
        out.boot.LOGODDS[i,1:NumGammaW] <- bs.results[[i]]$LOGODDS_SW
        out.boot.LOGODDS[i,NumGammaW+(1:NumGammaM)] <- bs.results[[i]]$LOGODDS_SM
      }
    }
#  bias-correct these
    coef.boot.bias <- apply(out.boot, 2, stats::median, na.rm=TRUE) - th_hat
#   recompute pmf_est
    # VIP using the uncorrected value to estimate pmf_est
    out$coefficients <- out$solution
    pmf_est <- exp(augpmfnew(out$coefficients[1:NumBeta],
                GammaW=out$coefficients[NumBeta+(1:NumGammaW)], 
                GammaM=out$coefficients[(NumBeta+NumGammaW)+(1:NumGammaM)],
                S, X, Z,
                pmfW, pmfM, gw=gw, gm=gm))
    pmf_est[nrow(pmf_est),ncol(pmf_est)] <- 0
    PMF_SW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
    out$PMF_SW <- PMF_SW / (PMF_SW + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],1,sum))
    PMF_SM <- pmf_est[ nrow(pmf_est),-ncol(pmf_est),drop=FALSE]
    out$PMF_SM <- PMF_SM / (PMF_SM + 0.5*apply(pmf_est[ -nrow(pmf_est),-ncol(pmf_est),drop=FALSE],2,sum))
    PMF_PW <- pmf_est[-nrow(pmf_est), ncol(pmf_est),drop=FALSE]
    out$PMF_PW <- 1 - out$PMF_SW
    out$PMF_PM <- 1 - out$PMF_SM 
    if(control$logodds_single){
      b <- sum(out$PMF_SW*pmfW)
      out$LOGODDS_SW <- out$coefficients[(NumBeta+1):(NumBeta+NumGammaW)] - log(b/(1-b))
      b <- sum(out$PMF_SM*pmfM)
      out$LOGODDS_SM <- out$coefficients[(NumBeta+NumGammaW+1):(NumBeta+NumGammaW+NumGammaM)] - log(b/(1-b))
    }else{
      out$LOGODDS_SW <- log(out$PMF_SW/(1-out$PMF_SW))
      out$LOGODDS_SM <- log(out$PMF_SM/(1-out$PMF_SM))
    }
#
    pmf_est[-nrow(pmf_est), -ncol(pmf_est)] <- 2*pmf_est[-nrow(pmf_est), -ncol(pmf_est)]
    pmf_est <- pmf_est/sum(pmf_est)
    out$coefficients <- th_hat-coef.boot.bias

    coef.boot.bias.LOGODDS <- apply(out.boot.LOGODDS, 2, stats::median, na.rm=TRUE) - th_hat[NumBeta+(1:NumGamma)]

    out$coefficients[NumBeta+(1:NumGammaW)] <- out$LOGODDS_SW
    out$coefficients[NumBeta+NumGammaW+(1:NumGammaM)] <- out$LOGODDS_SM

    out$ext.covar <- stats::cov(out.boot)
    if(control$robust.cov){
      pos.IQR <- apply(out.boot,2,stats::IQR,na.rm=TRUE)>0
      out.text <- capture.output(
        robust.cov <- try(MASS::cov.mcd(out.boot[,pos.IQR])),
       type="message")
      if(!any(startsWith(out.text,"Error in solve")) & !inherits(robust.cov,"try-error")){
        a <- cbind(sqrt(diag(out$ext.covar[pos.IQR,pos.IQR])),sqrt(diag(robust.cov$cov)))
        message_print(a)
        out$ext.covar[pos.IQR,pos.IQR] <- robust.cov$cov
      }
    }
    out$covar <- out$ext.covar[1:NumBeta,1:NumBeta]
	
    coef.boot <- out.boot[,1:NumBeta,drop=FALSE]
    theta_hat <- th_hat[1:NumBeta]
    theta_hat_SE <- apply(coef.boot, 2, stats::sd, na.rm=T)
    LB = control$alpha/2
    UB = 1-control$alpha/2
    theta_bootLB <- apply(coef.boot, 2, stats::quantile, probs = LB, na.rm=T)
    theta_bootUB <- apply(coef.boot, 2, stats::quantile, probs = UB, na.rm=T)
    out$basic.bootCI <- cbind(2*theta_hat-theta_bootUB,2*theta_hat-theta_bootLB)
    out$percent.bootCI <- cbind(theta_bootLB,theta_bootUB)
    out$theta_hat_SE <- theta_hat_SE
    out$t.bootCI <- cbind(theta_hat-1.96*theta_hat_SE, theta_hat+1.96*theta_hat_SE)
    a <- diag(out$ext.covar.hessian)[1:NumBeta]
    a[a<0] <- 0
    out$analyticalCI <- cbind(theta_hat-1.96*sqrt(a), theta_hat+1.96*sqrt(a))
    colnames(out$basic.bootCI) = colnames(out$percent.bootCI) = colnames(out$t.bootCI) = colnames(out$analyticalCI) = c('LB','UB')

    SE_boot <- sqrt(diag(out$covar))	
    tstar_LB <- apply(tstar.boot, 2, stats::quantile, probs = LB, na.rm=T)
    tstar_UB <- apply(tstar.boot, 2, stats::quantile, probs = UB, na.rm=T)
    out$student.bootCI <- c(theta_hat-tstar_UB*SE_boot, theta_hat-tstar_LB*SE_boot)
#  bias-correct these
    out.boot <- sweep(out.boot,2,2*coef.boot.bias,"-")
    out$boot.coef <- out.boot
    out.boot.LOGODDS <- sweep(out.boot.LOGODDS,2,2*coef.boot.bias.LOGODDS,"-")
    out$boot.coef.LOGODDS <- out.boot.LOGODDS
    }else{ # non-bootstrap version
      out$coefficients <- out$solution
    } # end of the bootstrap

    if(!is.na(match("Int", colnames(Xdata)))){
      out$Xdata <- Xdata[, colnames(Xdata) != "Int"]
    }else{
      out$Xdata <- Xdata
    }
    if(!is.na(match("Int", colnames(Zdata)))){
      out$Zdata <- Zdata[, colnames(Zdata) != "Int"]
    }else{
      out$Zdata <- Zdata
    }
    out$NumBeta <- NumBeta
    out$NumBetaW <- NumBetaW
    out$NumBetaM <- NumBetaM
    out$NumGamma  <- NumGamma
    out$NumGammaW <- NumGammaW
    out$NumGammaM <- NumGammaM
    out$num_women <- num_women
    out$num_men <- num_men
    out$nobs <- num_sampled
    
    out$pmf <- pmf
    out$pmfW=pmfW; out$pmfM=pmfM
    out$pmfN=pmfN
    out$pmf_est <- pmf_est
    
    # For GOF
    out$bs.results <- bs.results

    class(out) <- "rpm"
    
    return(out)
}
