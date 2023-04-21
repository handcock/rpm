#' Summarize Revealed Preference Matchings data via a Model Specification
#' 
#' \code{\link{summary_rpm}} produces tabular summaries of data revealed preference matchings
#' based on a formula specifying a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' 
#' @aliases print.summary_rpm show.summary_rpm
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
#' \item{X_w}{: This is character string of the name of the weight variable for women. The sum of the weights should be the
#' number of women in the population.}
#' \item{Z_w}{: This is character string of the name of the weight variable for men. The sum of the weights should be the
#' number of men in the population.}
#' \item{pair_w}{: This is character string of the name of the weight variable for pairs.}
#' }
#' @return \code{\link{summary}} returns a list with many components, like \code{\link{rpm}} object without the model estimates. In particular it includes \code{stats} and \code{popstats}.
#' \code{stats} is the named vector of sample statistics from the model. 
#' while \code{popstats} is the named vector of population statistics from the model. 
#' It alos includes \code{counts} and \code{pmf}. Each of these is a contingency table in array
#' representation of S3 class \code{c("xtabs", "table")}, with a \code{"call"}
#' @seealso control.rpm, summary.rpm, rpm
#' @references
#'
#' Goyal, Shuchi; Handcock, Mark S.; Jackson, Heide M.; Rendall, Michael S. and Yeung, Fiona C. (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#'
#' Dagsvik, John K. (2000) \emph{Aggregation in Matching Markets} \emph{International Economic Review}, Vol. 41, 27-57.
#' JSTOR: https://www.jstor.org/stable/2648822, \doi{10.1111/1468-2354.00054}
#'
#' Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941. \doi{10.3982/ECTA12299}
#' @keywords models
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' summary_rpm(~match("edu") + WtoM_diff("edu",3),
#'         Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'         X_w="X_w", Z_w="Z_w",
#'         pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'         sampled="sampled",sampling_design="stock-flow")
#' @export summary_rpm
summary_rpm <- function(formula, Xdata, Zdata,
    Xid=NULL, Zid=NULL, pair_id=NULL,
    X_w=NULL, Z_w=NULL, pair_w=NULL,
    sampled=NULL,
    sampling_design="stock-flow",
    control=control.rpm(), verbose=FALSE){

    if(!is.null(control$seed)) set.seed(control$seed)

    if(missing(Xid)){
      warning("The variable name in Xdata of the women's IDs is missing and is presumed to be 'Xid'.")
      Xid <- "Xid"
      if(is.null(Xdata[["Xid"]])){
        stop("The variable name in Xdata of the women's IDs must be specified as the character string 'Xid'.")
      }
    }
    if(missing(Zid)){
      warning("The variable name in Zdata of the men's IDs is missing and is presumed to be 'Zid'.")
      Zid <- "Zid"
      if(is.null(Zdata[["Zid"]])){
        stop("The variable name in Zdata of the women's IDs must be specified as the character string 'Xid'.")
      }
    }
    if(missing(pair_id)){
      warning("The variable name in Xdata of the men's IDs of the women's partners has not been specified. It is presumed to be 'pair_id'.")
      pair_id <- "pair_id"
      if(is.null(Xdata[["pair_id"]])){
        stop("The variable 'pair_id' in Xdata of the men's IDs of the women's partners is NULL.")
      }
    }
    if(missing(sampled) & sampling_design != "census"){
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
      print(a)
      message("Please correct this before rerunning.")
    }

    if(!all(M_paired_to_W)){
      message("The men's data, Zdata, contains people with pair IDs that do not appear in the women's data, Xdata. These are:")
      a <- Zdata[!is.na(Zdata[,pair_id]),][!M_paired_to_W,]
      print(a)
      message("Please correct this before rerunning.")
    }

    if(!is.empty(W_duplicated_pair_id)){
      message("The women's data, Xdata, contains duplicated pair IDs. These are:")
      a <- Xdata[!is.na(Xdata[,pair_id]),][W_pair_id %in% W_duplicated_pair_id,]
      print(a[order(a[,pair_id]),])
      message("Please correct this before reruning.")
    }
    if(!is.empty(M_duplicated_pair_id)){
      message("The men's data, Zdata, contains duplicated pair IDs. These are:")
      a <- Zdata[!is.na(Zdata[,pair_id]),][M_pair_id %in% M_duplicated_pair_id,]
      print(a[order(a[,pair_id]),])
      message("Please correct this before reruning.")
    }
    if(!is.empty(W_duplicated_pair_id) | !is.empty(W_duplicated_pair_id) | !all(W_paired_to_M) | !all(M_paired_to_W)){
      stop()
    }

    if(sampling_design != "census"){
    # IDs of the women matched to the sampled men (and vice versa)
     paired_and_sampled_M <- Zdata[,sampled] & !is.na(Zdata[,pair_id])
     paired_and_sampled_W <- Xdata[,sampled] & !is.na(Xdata[,pair_id])
     paired_and_unsampled_M <- !Zdata[,sampled] & !is.na(Zdata[,pair_id])
     paired_and_unsampled_W <- !Xdata[,sampled] & !is.na(Xdata[,pair_id])
     W_paired_to_sampled_M <- match(Zdata[paired_and_sampled_M,pair_id], Xdata[,Xid])
     M_paired_to_sampled_W <- match(Xdata[paired_and_sampled_W,pair_id], Zdata[,Zid])
     W_paired_to_unsampled_M <- match(Zdata[paired_and_unsampled_M,pair_id], Xdata[,Xid])
     M_paired_to_unsampled_W <- match(Xdata[paired_and_unsampled_W,pair_id], Zdata[,Zid])
    }else{
     paired_and_sampled_M <- !is.na(Zdata[,pair_id])
     paired_and_sampled_W <- !is.na(Xdata[,pair_id])
     W_paired_to_sampled_M <- match(Zdata[paired_and_sampled_M,pair_id], Xdata[,Xid])
     M_paired_to_sampled_W <- match(Xdata[paired_and_sampled_W,pair_id], Zdata[,Zid])
    }

    if(is.null(pair_w) & sampling_design != "census"){
    # This is the unweighted case
    # Construct individual weight case 
#     a=rep(NA, length=nrow(Xdata))
      # sampled women matched to men: women + man
#     a[paired_and_sampled_W] <- 1/(1/Xdata[paired_and_sampled_W,X_w] + 1/Zdata[M_paired_to_sampled_W,Z_w])
      # unsampled women matched to men: women + man
#     a[!Xdata[,sampled] & !is.na(Xdata[,pair_id])] <- 1/(1/Zdata[match(Xdata[!Xdata[,sampled] & !is.na(Xdata[,pair_id]),Xid],Zdata[,pair_id]),Z_w] + 1/Xdata[W_paired_to_sampled_M,X_w])
#     a[paired_and_unsampled_W] <- 1/(1/Zdata[match(Xdata[paired_and_unsampled_W,Xid],Zdata[,pair_id]),Z_w] +
#                                     1/Xdata[paired_and_unsampled_W,X_w])
      a=rep(NA, length=nrow(Xdata))
      a[paired_and_sampled_W] <- Xdata[paired_and_sampled_W,X_w]
#     a[paired_and_unsampled_W] <- Zdata[match(Xdata[paired_and_unsampled_W,Xid],Zdata[,pair_id]),Z_w]
      a[paired_and_unsampled_W] <- 0
#     a[paired_and_unsampled_W] <- Xdata[paired_and_unsampled_W,X_w]
      pair_w <- "pair_w"
      Xdata[[pair_w]]=a
#     a=rep(NA, length=nrow(Zdata))
#     a[paired_and_sampled_M] <- 1/(1/Zdata[paired_and_sampled_M,Z_w] + 1/Xdata[W_paired_to_sampled_M,X_w])
#     a[paired_and_unsampled_M] <- 1/(1/Xdata[match(Zdata[paired_and_unsampled_M,Zid],Xdata[,pair_id]),X_w] +
#                                     1/Zdata[paired_and_unsampled_M,Z_w])
      a=rep(NA, length=nrow(Zdata))
      a[paired_and_sampled_M] <- Zdata[paired_and_sampled_M,Z_w]
#     a[paired_and_unsampled_M] <- Xdata[match(Zdata[paired_and_unsampled_M,Zid],Xdata[,pair_id]),X_w]
      a[paired_and_unsampled_M] <- 0
      Zdata[[pair_w]]=a
    }
    if(is.null(pair_w)){ pair_w <- "pair_w" }
    
    # get the proportion of men and women
    
    if(sampling_design != "census"){
     if (sampling_design == "stock-stock") {
      XdataM <- Xdata[W_paired_to_sampled_M,]
      ZdataW <- Zdata[M_paired_to_sampled_W,]
      # Presumes individual weights
      n_w = sum(Xdata[Xdata[,sampled] & is.na(Xdata[,pair_id]),X_w])
      n_w = n_w + sum(Xdata[Xdata[,sampled] & !is.na(Xdata[,pair_id]),X_w]) # number of women
      n_m = sum(Zdata[Zdata[,sampled] & is.na(Zdata[,pair_id]),Z_w])
      n_m = n_m + sum(Zdata[Zdata[,sampled] & !is.na(Zdata[,pair_id]),Z_w]) # number of men
      n = n_w + n_m
      gw = log(n_w/n) # to ensure exp(gw)+exp(gm) = 1
      gm = log(n_m/n) # to ensure exp(gw)+exp(gm) = 1
     }else{
      # Presumes individual weights
      n_w = sum(Xdata[Xdata[,sampled],X_w])
      n_m = sum(Zdata[Zdata[,sampled],Z_w])
      n = n_w + n_m
      gw = log(n_w/n) # to ensure exp(gw)+exp(gm) = 1
      gm = log(n_m/n) # to ensure exp(gw)+exp(gm) = 1
     }
    }else{
      n = nrow(Xdata) + nrow(Zdata) # The population size
      gw = log(nrow(Xdata)/n) # to ensure exp(gw)+exp(gm) = 1
      gm = log(nrow(Zdata)/n)
    }
 
    # 1) parse the formula
    # intercept is always added as the first column
    Xdata <- cbind(1, Xdata)
    Zdata <- cbind(1, Zdata)
    colnames(Xdata)[1] <- "Int"
    colnames(Zdata)[1] <- "Int"
    
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
      print(model_vars[!W_matched_vars])
      message("Please correct this before reruning.")
    }
    M_matched_vars <- model_vars %in% colnames(Zdata)    
    if(!all(M_matched_vars)){
      message("Some of the variables in the formula do not appear in the men's data, Zdata. These are:")
      print(model_vars[!M_matched_vars])
      message("Please correct this before reruning.")
    }
    if(!all(W_matched_vars) | !all(W_matched_vars)){
      stop()
    }

    # 3) Compute the marginal distributions of women's and men's types
    Xu <- unique(Xdata[,model_vars])
    Xu <- Xu[do.call(order, as.data.frame(Xu)),]
    Zu <- unique(Zdata[,model_vars])
    Zu <- Zu[do.call(order, as.data.frame(Zu)),]

    cnW <- paste(colnames(Xu)[2],Xu[,2], sep=".")
    if(ncol(Xu)>2){
    for(i in 3:ncol(Xu)){
      cnW <- paste(cnW,paste(colnames(Xu)[i],Xu[,i], sep="."),sep='.')
    }}
    cnM <- paste(colnames(Zu)[2],Zu[,2], sep=".")
    if(ncol(Zu)>2){
    for(i in 2:ncol(Zu)){
      cnM <- paste(cnM,paste(colnames(Zu)[i],Zu[,i], sep="."),sep='.')
    }}
    
    # 4) Create joint PMF
    # Xtype: group membership for women (one for each woman in the pop)
    Xtype <- rep(NA,nrow(Xdata))
    for(i in 1:nrow(Xu)){
        Xtype[apply(Xdata[,model_vars], 1, function(x) identical(x, unlist(Xu[i,])))] <- i
    }
    Xdata[["Xtype"]] <- Xtype
    rm(Xtype)
    # Ztype: group membership for men (one for each man in the pop)
    Ztype <- rep(NA,nrow(Zdata))
    for(i in 1:nrow(Zu)){
        Ztype[apply(Zdata[,model_vars], 1, function(x) identical(x, unlist(Zu[i,])))] <- i
    }
    Zdata[["Ztype"]] <- Ztype
    rm(Ztype)
    
    if (!is.null(X_w) & sampling_design != "census") {
        Xdata$X_w = Xdata[,X_w]
    }
    if (!is.null(Z_w) & sampling_design != "census") {
        Zdata$Z_w = Zdata[,Z_w]
    }

    num_Xu = nrow(Xu)
    num_Zu = nrow(Zu)
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
    # Compute the number sampled (manly for s.e. computation
    if (sampling_design == "stock-stock") {
      num_sampled = sum(Xdata[,sampled] & is.na(Xdata[,pair_id]))
      num_sampled = num_sampled + sum(Zdata[,sampled])
    }else{
      if (sampling_design == "stock-flow") {
        num_sampled <- sum(Xdata[,sampled])+sum(Zdata[,sampled])
      }else{
        num_sampled <- nrow(Xdata) + nrow(Zdata)
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

      # The number of pairs of people in the sample of each (X,Z) pair
      counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      
      if (!is.empty(Xcounts_single)) {
        pmf[1:num_Xu,1+num_Zu] = Xcounts_single / N
        counts[1:num_Xu,1+num_Zu] = Xcounts_single
      }
      if (!is.empty(Zcounts_single)) {
        pmf[1+num_Xu,1:num_Zu] = Zcounts_single / N
        counts[1+num_Xu,1:num_Zu] = Zcounts_single
      }
      
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
      pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      # The proportion of people in the population of each (X,Z) pair
      pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N

      # The number of pairs of people in the sample of each (X,Z) pair
      counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu)))
      
      if (!is.empty(Xtype_single)) {
        pmf[1:num_Xu,1+num_Zu] = Xtype_single / N
        counts[1:num_Xu,1+num_Zu] = Xcounts_single
      }
      if (!is.empty(Ztype_single)) {
        pmf[1+num_Xu,1:num_Zu] = Ztype_single / N
        counts[1+num_Xu,1:num_Zu] = Zcounts_single
      }

    } else { # assume "stock-flow"
      
      pmf = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      counts = matrix(0,nrow=1+num_Xu, ncol=1+num_Zu) # women (X) indexed by row, men (Z) indexed by column
      colnames(pmf) <- c(cnM,"singles")
      rownames(pmf) <- c(cnW,"singles")
      colnames(counts) <- c(cnM,"singles")
      rownames(counts) <- c(cnW,"singles")
      
      # The number of people in the population
      N <- n

      # The number of people in the population of each (X,Z) pair
      pmf[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(Xdata[paired_and_sampled_W,X_w]~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(Zdata[paired_and_sampled_M,Z_w]~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      # The proportion of people in the population of each (X,Z) pair
      # note that pmf[1:num_Xu,1:num_Zu] is *twice* the number of pairs so that
      # exp(-gw)*apply(pmf,1,sum) = pmfW
      pmf[1:num_Xu,1:num_Zu] <- pmf[1:num_Xu,1:num_Zu] / N

      # The number of pairs of people in the sample of each (X,Z) pair
      counts[1:num_Xu,1:num_Zu] <- as.numeric(stats::xtabs(~factor(Xdata[paired_and_sampled_W,"Xtype"],1:num_Xu)+factor(ZdataW[,"Ztype"],1:num_Zu))) + as.numeric(stats::xtabs(~factor(XdataM[,"Xtype"],1:num_Xu)+factor(Zdata[paired_and_sampled_M,"Ztype"],1:num_Zu)))
      
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

    }

    pmfN <- pmf*num_sampled

    if(verbose){
      message(sprintf("Population proportions of women by category:"))
      print(pmfW)
      message(sprintf("Population proportions of men by category:"))
      print(pmfM)
      message(sprintf("Matrix of population proportions of women x men by category:"))
      print(pmf)
    }

    out <- list(control=control)

    modelmat <- rpm.model.matrix(model.terms, Xu, Zu)
 
    if(length(modelmat$S)>0){
      Sstats = as.vector(apply(modelmat$S,3,function(x){sum(x*counts[-nrow(counts),-ncol(counts)])}))
      names(Sstats) <- modelmat$Snames
      Spopstats = as.vector(apply(modelmat$S,3,function(x){sum(x*pmf[-nrow(pmf),-ncol(pmf)])}))
      names(Spopstats) <- modelmat$Snames
    }else{
      Sstats <- NULL
      Spopstats <- NULL
    }
    if(length(modelmat$X)>0){
      Xstats = as.vector(apply(modelmat$X,3,function(x){sum(x*counts[-nrow(counts),-ncol(counts)])}))
      names(Xstats) <- modelmat$Xnames
      Xpopstats = as.vector(apply(modelmat$X,3,function(x){sum(x*pmf[-nrow(pmf),-ncol(pmf)])}))
      names(Xpopstats) <- modelmat$Xnames
    }else{
      Xstats <- NULL
      Xpopstats <- NULL
    }
    if(length(modelmat$Z)>0){
      Zstats = as.vector(apply(modelmat$Z,c(1,2),function(x){sum(x*counts[-nrow(counts),-ncol(counts)])}))
      names(Zstats) <- modelmat$Znames
      Zpopstats = as.vector(apply(modelmat$Z,c(1,2),function(x){sum(x*pmf[-nrow(pmf),-ncol(pmf)])}))
      names(Zpopstats) <- modelmat$Znames
    }else{
      Zstats <- NULL
      Zpopstats <- NULL
    }
 
    out$stats <- c(Sstats, Xstats, Zstats)
    out$popstats <- c(Spopstats, Xpopstats, Zpopstats)

    class(counts) <- c("xtabs", "table")
    attr(counts, "call") <- match.call()

    out$pmf=pmf; out$pmfW=pmfW; out$pmfM=pmfM
    out$counts=counts
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
    out$Xdata <- Xdata
    out$Zdata <- Zdata
    out$gw <- gw
    out$gm <- gm
    out$N <- N
    out$sampling_design <- sampling_design
    
    class(out) <- "summary_rpm"
    
    return(out)
}
#' @method print summary_rpm
#' @export
# The <print.summary_rpm> function prints a subset of the information given
# by <summary_rpm>
print.summary_rpm <- function (x, 
                              digits = max(6, getOption("digits") - 3),
                              ...){
  
    cat("\n==========================\n")
    cat("Summary of RPM model data\n")
    cat("==========================\n\n")
  
    cat("Formula:   ")
    print(x$formula)
    cat("\n")
  
    cat("\n==========================\n")
    cat("Sample statistics\n")
    cat("==========================\n\n")
    print(x$stats)

    cat("\n==========================\n")
    cat("Population statistics\n")
    cat("==========================\n\n")
    print(x$popstats)

    message(sprintf("\nSampling Design: %s",x$sampling_design))
    message(sprintf("Sample size: %f",x$num_sampled))
    message(sprintf("Matrix of counts:"))
    
    print(x$counts)

    message(sprintf("\nPopulation size: %f",x$N))
    message(sprintf("\nPopulation proportions of women by category:"))
    print(x$pmfW)
    message(sprintf("\nPopulation proportions of men by category:"))
    print(x$pmfM)
    message(sprintf("\nMatrix of population proportions of women x men by category:"))
    print(round(x$pmf,digits))

  invisible(x)
}

#' @method print summary.rpm
#' @export
show.summary_rpm <- print.summary_rpm
