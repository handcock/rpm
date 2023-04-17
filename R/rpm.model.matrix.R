#' Creates a model matrix to estimate the parameters of a Revealed Preference Matchings Model
#' 
#' \code{\link{rpm.model.matrix}} assumes a bipartite network (i.e. two-sided matching market)
#' It creates a model matrix according to the formula passed in.
#' See \code{\link{rpm-terms}} for a description of the possible terms.
#' 
#' @param model.terms For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}. This includes the 
#' covariates used to construct the model matrix.
#' They are used in conjunction with the model terms. 
#' @param Xall the unique types of women
#' @param Zall the unique types of men
#' @param intercept logical; If TRUE, the default, an intercept term is prepended.
#' @return A list consists of the following elements:
#' \item{X}{the model matrix for women.}
#' \item{Z}{the model matrix for men.}
#' \item{Xnames}{the names of the covariates for women.} 
#' \item{Znames}{the names of the covariates for men.}
#' @seealso rpm
#' @references Goyal, Handcock, Jackson. Rendall and Yeung (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#' Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' # nothing yet
#' @export rpm.model.matrix
#' 
rpm.model.matrix = function(model.terms, Xall, Zall, intercept=TRUE)
{
    # assumes Xall and Zall are the unique types
    
    # no. of members on each side
    nW = nrow(Xall) 
    nM = nrow(Zall)

    ncov = length(model.terms)
    
    # returned variables
    S = array(0, dim=c(nW,nM,0))
    X = array(0, dim=c(nW,nM,0))
    Z = array(0, dim=c(nM,nW,0))
    Snames = NULL
    Xnames = NULL
    Znames = NULL
    
    for(i in 1:ncov)
    {
        switch(as.character(model.terms[[i]][[1]]),
               '1' = {
                   if(attr(model.terms,"sign")[[i]]>0){
                     S = abind::abind(S,
                           array(1, dim=c(nW,nM,1)),
                           along = 3)
                     Snames = c(Snames, "intercept")
                   }
                   intercept <- FALSE
               },
               absdiff = {
                   
                   attrname = as.character(model.terms[[i]][[2]])
                   if(length(model.terms[[i]])<3){
                     S = abind::abind(S,
                           1*outer(Xall[,attrname], Zall[,attrname], function(x,z){1*(abs(x-z))}),  
                           along = 3)
                     Snames = c(Snames, "absdiff")
                   }else{
                     difflev = as.numeric(model.terms[[i]][[3]])
                     S = abind::abind(S,
                           1*outer(Xall[,attrname], Zall[,attrname], function(x,z){1*(abs(x-z)==difflev)}),  
                           along = 3)
                     Snames = c(Snames, paste("absdiff", attrname, difflev, sep="."))
                   }
               },
               W_cov = {
                   for(j in 2 : length(model.terms[[i]]))
                   {	
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       # every row has the same value -- group 1's raw values for this attribute
                       X = abind::abind(X,
                                 matrix(rep(Xall[, attrname], times = nM), ncol = nM, byrow = FALSE),
                                 along = 3)
                       Xnames = c(Xnames, paste("W_cov", attrname, sep="."))
                   }
               },
               M_cov = {
                   for(j in 2 : length(model.terms[[i]]))
                   {	
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       # every row has the same value -- group 2's raw values for this attribute
                       Z = abind::abind(Z,
                                 matrix(rep(Zall[, attrname], times = nW), ncol = nW, byrow = FALSE),
                                 along = 3)
                       Znames = c(Znames, paste("M_cov", attrname, sep="."))
                   }
               },
               W_factor = {
                   for(j in 2 : length(model.terms[[i]]))
                   {	
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       # get all the factor levels possible for both sides
                       temp = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       for(k in temp[-1])
                       {
                           # every row has the same value --  1 if group 1's value for this attribute is same as k; 0 otherwise
                           X = abind::abind(X, 
                                     matrix(rep((1 * (Xall[, attrname] == k)), times = nM), ncol = nM, byrow = FALSE),
                                     along = 3)
                           Xnames = c(Xnames, paste("W_factor", attrname, k, sep="."))
                       }
                   }
               },
               M_factor = {
                   for(j in 2 : length(model.terms[[i]]))
                   {			
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       # get all the factor levels possible for both sides
                       temp = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       for(k in temp[-1])
                       {
                           # every row has the same value --  1 if group 2's value for this attribute is same as k; 0 otherwise
                           Z = abind::abind(Z, 
                                     matrix(rep((1 * (Zall[, attrname] == k)), times = nW), ncol = nW, byrow = FALSE),
                                     along = 3)
                           Znames = c(Znames, paste("M_factor", attrname, k, sep="."))
                       }
                   }
               },
               homophily = {
                   parts = as.character(model.terms[[i]][[2]])
                   
                   # only the matches for the main effect 
                   if(length(parts) == 1)
                   {
                       attrname = parts
                    
                       # matching values on attribute
                       S = abind::abind(S,
                                  1*outer(Xall[, attrname], Zall[, attrname], "=="),
                                  along = 3) 
                       Snames = c(Snames, paste("homophily", attrname, sep="."))
                   }
                   else # matching on more than a single attribute
                   {
                       unique.combo = unique(rbind(Xall[,parts], Zall[,parts]))
                       
                       # assign unique group membership for group 1 
                       Xgroupmem = rep(NA,nrow(Xall))
                       for(i in 1:nrow(unique.combo)){
                           Xgroupmem[apply(Xall[,parts], 1, function(x) identical(x, unique.combo[i,]))] = i
                       }
                       
                       # assign unique group membership for group 2 
                       Zgroupmem = rep(NA,nrow(Zall))
                       for(i in 1:nrow(unique.combo)){
                           Zgroupmem[apply(Zall[,parts], 1, function(x) identical(x, unique.combo[i,]))] = i
                       }
                       
                       # matching values on attribute
                       S = abind::abind(S,
                                  1*outer(Xgroupmem, Zgroupmem, "=="),
                                  along = 3) 
                       Snames = c(Snames, paste("Whomophily", paste(parts, sep=".", collapse = "."), sep="."))
                   }
               },
               match = {
                   attrname = as.character(model.terms[[i]][[2]])
                   
                   # only the matches for the main effect 
                   if("diff" %in% names(model.terms[[i]]))
                   {
                       # get all the factor levels possible for both sides
                       u = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       if("collapse" %in% names(model.terms[[i]])){
                         cl <- eval(model.terms[[i]][["collapse"]],envir=list())
                         for(k in seq_along(cl)){
                           del <- rep(FALSE,length(u))
                           for (m in 1:length(u)){
                            del[m] = (u[m] %in% cl[[k]])
                           }
                           u = u[!del]
                         }
                       }else{
                         cl <- NULL
                       }
                       if(length(u)>0){
                        cll <- NULL
                        for (m in 1:length(u)){
                         cll <- c(cll,list(u[m]))
                        }
                        cl <- c(cll,cl)
                       }

                       # matching on each level
                       for (k in seq_along(cl))
                       {
                           # only matching on kth level
                           W = rep(FALSE, length(Xall[,attrname]))
                           for (l in seq_along(cl[[k]])){
                              W = W | ((Xall[,attrname] == cl[[k]][[l]]) & (Zall[,attrname] == cl[[k]][[l]]))
                           }
                           
                           S = abind::abind(S,
                                      1*diag(W),
                                      along = 3) 
                           Snames = c(Snames, paste("match", attrname, paste0(cl[[k]],collapse="."), sep=".",collapse=""))
                       }
                   }
                   else # uniform matching 
                   {
                       # get all the factor levels possible for both sides
                       u = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       # matching on each level
                       for (k in u)
                       {   
                           # only matching on kth level
                           W = 1 * (Xall[,attrname] == k)
                           M = 1 * (Zall[,attrname] == k)
                           
                           S = abind::abind(S,
                                      1*outer(W, M, "*"),
                                      along = 3) 
                           Snames = c(Snames, paste("match", attrname, k, sep="."))
                       }
                   }
               },
               mix = {
                   attrname = as.character(model.terms[[i]][[2]])
                   # get all the factor levels possible for both sides
                   b1 = sort(unique(Xall[,attrname]))
                   b2 = sort(unique(Zall[,attrname]))
                   unique.combo <- expand.grid(row = b1, col = b2, stringsAsFactors=FALSE)
                       
                   if("base" %in% names(model.terms[[i]])){
                     base <- eval(model.terms[[i]][["base"]],envir=list())
                     if(is.null(base)){
                       sel <- 1:nrow(unique.combo)
                     }else{
                       sel <- (1:nrow(unique.combo))[-base]
                     }
                   }else{
                     sel <- 2:nrow(unique.combo)
                   }
                   bcombo <- unique.combo[sel,]
                   
                   if("collapse" %in% names(model.terms[[i]])){
                     cl <- eval(model.terms[[i]][["collapse"]],envir=list())
                     for(k in seq_along(cl)){
                      for (l in seq_along(cl[[k]])){
                       del <- rep(FALSE,nrow(bcombo))
                       for (m in 1:nrow(bcombo)){
                        del[m] = (bcombo[m,1] == cl[[k]][[l]][1]) & (bcombo[m,2] == cl[[k]][[l]][2])
                       }
                       bcombo = bcombo[!del,]
                     }
                    }
                   }else{
                     cl <- NULL
                   }
                   cll <- NULL
                   for (m in 1:nrow(bcombo)){
                     cll <- c(cll, list(list(bcombo[m,])))
                   }
                   cl <- c(cll, cl)
                    
                   # matching on each combination of level
                   for (j in seq_along(cl))
                   {
                        Sn <- paste("mix", attrname, sep='.')
                        Sj <- matrix(0,ncol=length(Zall[,attrname]), nrow=length(Xall[,attrname]))
                        for (l in seq_along(cl[[j]])){
                           Sj[Xall[,attrname] %in% cl[[j]][[l]][1], Zall[,attrname] %in% cl[[j]][[l]][2]] <- 1
                           Sn <- paste(Sn, paste0(cl[[j]][[l]],collapse="."),sep='_')
                        }

                        S = abind::abind(S,
                                   Sj,
                                   along = 3) 
                        Snames = c(Snames, Sn)
                   }
               },
               W_smallerthan = {
                   for(j in 2 : length(model.terms[[i]]))
                   {			
                       attrname = as.character(model.terms[[i]][[j]])
                   
                       X = abind::abind(X,
                                 1*outer(Xall[,attrname], Zall[,attrname],">"),  # a liking for group 2 of lower level
                                 along = 3)
                       Xnames = c(Xnames, paste("W_smallerthan", attrname, sep="."))
                   }
               },
               W_greaterthan = {
                   for(j in 2 : length(model.terms[[i]]))
                   {			
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       X = abind::abind(X,
                                 1*outer(Xall[,attrname], Zall[,attrname],"<"),  # a liking for group 2 of higher level
                                 along = 3)
                       Xnames = c(Xnames, paste("W_greaterthan", attrname, sep="."))
                   }
               },
               M_smallerthan = {
                   for(j in 2 : length(model.terms[[i]]))
                   {			
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       Z = abind::abind(Z,
                                 1*outer(Zall[,attrname], Xall[,attrname],">"),  # a liking for group 1 of lower level
                                 along = 3)
                       Znames = c(Znames, paste("M_smallerthan", attrname, sep="."))
                   }
               },
               M_greaterthan = {
                   for(j in 2 : length(model.terms[[i]]))
                   {			
                       attrname = as.character(model.terms[[i]][[j]])
                       
                       Z = abind::abind(Z,
                                 1*outer(Zall[,attrname], Xall[,attrname],"<"),  # a liking for group 1 of higher level
                                 along = 3)
                       Znames = c(Znames, paste("M_greaterthan", attrname, sep="."))
                   }
               },
               diff = {
                   
                   if(length(model.terms[[i]])<2){
                     stop("The 'diff' term requires a first argument, the attribute to take the difference of.")
                   } 
                   attrname = as.character(model.terms[[i]][[2]])
                       
                   X = abind::abind(X,
                             1*outer(Xall[,attrname], Zall[,attrname], function(x,z){x-z}),  
                             along = 3)
                   Xnames = c(Xnames, paste("diff", attrname, sep="."))
               },
               WtoM_diff = {
                   
                   if(length(model.terms[[i]])<2){
                     stop("The 'WtoM_diff' term requires a first argument, the attribute to take the difference of.")
                   } 
                   attrname = as.character(model.terms[[i]][[2]])
                   if(length(model.terms[[i]])<3){
                     stop("The 'WtoM_diff' term requires a second argument, the difference level.")
                   } 
                   difflev = as.numeric(model.terms[[i]][[3]])
                       
                   X = abind::abind(X,
                             1*outer(Xall[,attrname], Zall[,attrname], function(x,z){1*((x-z)==difflev)}),  
                             along = 3)
                   Xnames = c(Xnames, paste("WtoM_diff", attrname, difflev, sep="."))
               },
               MtoW_diff = {
                   
                   if(length(model.terms[[i]])<2){
                     stop("The 'MtoW_diff' term requires a first argument, the attribute to take the difference of.")
                   } 
                   attrname = as.character(model.terms[[i]][[2]])
                   if(length(model.terms[[i]])<3){
                     stop("The 'WtoM_diff' term requires a second argument, the difference level.")
                   } 
                   difflev = as.numeric(model.terms[[i]][[3]])
                   
                   Z = abind::abind(Z,
                             1*outer(Zall[,attrname], Xall[,attrname], function(z,x){1*((z-x)==difflev)}),  
                             along = 3)
                   Znames = c(Znames, paste("MtoW_diff", attrname, difflev, sep="."))
               },
               women = {
                   w = array(0, dim=c(nW,nM))
                   for(i in 1:nW){
                     w[i,] <- 1
                     X = abind::abind(X, w, along = 3)
                     w[i,] <- 0
                   }
                   Xnames = c(Xnames, paste("woman", 1:nW, sep="."))
               },
               men = {
                   m = array(0, dim=c(nW,nM,nM))
                   for(i in 1:nM){
                     m[,i,i] <- 1
                   }
                   X = abind::abind(X, m, along = 3)
                   Xnames = c(Xnames, paste("man", 1:nM, sep="."))
               },
               # default
               {
                   message("Error: unrecognized rpm term: ", as.character(model.terms[[i]][2]))
               }
               
               )
    }
    if(intercept){
      S = abind::abind(array(1, dim=c(nW,nM,1)), S, along = 3)
      Snames = c("intercept", Snames)
    }

    return(list(S=S ,X=X, Z=Z, Snames=Snames, Xnames=Xnames, Znames=Znames))
}
