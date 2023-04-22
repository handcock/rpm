#' Creates a model function list for the continuous terms in a Revealed Preference Matchings Model
#' 
#' \code{\link{rpm.model.matrix}} assumes a bipartite network (i.e. two-sided matching market)
#' It creates a model matrix according to the formula passed in.
#' See \code{\link{rpm-terms}} for a description of the possible terms.
#' 
#' @param model.terms For the details on the possible countinuous \code{<model terms>}, see
#' \code{\link{rpm-terms}}. This includes the 
#' covariates used to construct the model matrix.
#' They are used in conjunction with the model terms. 
#' @param control A list of control parameters for algorithm tuning. Constructed using
#' \code{\link{control.rpm}}, which should be consulted for specifics. 
#' @return A list of model terms as bivariate functions.
#' @seealso rpm
#' @references
#'
#' Goyal, Shuchi; Handcock, Mark S.; Jackson, Heide M.; Rendall, Michael S. and Yeung, Fiona C. (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation},
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 
#'
#' Dagsvik, John K. (2000) \emph{Aggregation in Matching Markets} \emph{International Economic Review},, Vol. 41, 27-57.
#' JSTOR: https://www.jstor.org/stable/2648822, \doi{10.1111/1468-2354.00054}
#'
#' Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941. \doi{10.3982/ECTA12299}
#' @keywords models
#' @examples
#' # nothing yet
#' @export rpm.model.functions
#' 
rpm.model.functions = function(model.terms, control)
{
    # assumes Xall and Zall are the unique types
    
    ncov = length(model.terms)
    
    # returned variables
    S = Sr = Spaired = list()
    Snames = AlphaS_K = NULL
    
    for(i in 1:ncov)
    {
        switch(as.character(model.terms[[i]][[1]]),
               absdiff = {
                   
                   S = c(S, 
                         function(x,z,K=NULL,xr=NULL,zr=NULL){
                          x = as.matrix(x); z = as.matrix(z); 
                          if(!is.null(xr)){
                            x = sweep(sweep(x,2,xr[2,]-xr[1,],"*"),2,xr[1,],"+")
                            z = sweep(sweep(z,2,zr[2,]-zr[1,],"*"),2,zr[1,],"+")
                          }
                          as.vector(1*(abs(x-z)))
                         }
                        )
                   Sr = c(Sr, 
                         function(x,z,K=NULL,xr=NULL,zr=NULL){
                          x = as.matrix(x); z = as.matrix(z); 
                          if(!is.null(xr)){
                            x = sweep(sweep(x,2,xr[2,]-xr[1,],"*"),2,xr[1,],"+")
                            z = sweep(sweep(z,2,zr[2,]-zr[1,],"*"),2,zr[1,],"+")
                          }
                          as.vector(1*(abs(x-z)))
                         }
                        )
                   Spaired = c(Spaired, 
                         function(x,z,K=NULL,xr=NULL,zr=NULL){
                          x = as.matrix(x); z = as.matrix(z); 
                          if(!is.null(xr)){
                            x = sweep(sweep(x,2,xr[2,]-xr[1,],"*"),2,xr[1,],"+")
                            z = sweep(sweep(z,2,zr[2,]-zr[1,],"*"),2,zr[1,],"+")
                          }
                          as.vector(1*(abs(x-z)))
                         }
                        )
                   AlphaS_K = c(AlphaS_K, 1)
                   Snames = c(Snames, "absdiff")
               },
               WtoM_diff = {
                   
                   if(length(model.terms[[i]])<2){
                     stop("The 'WtoM_diff' term requires a first argument, the attribute to take the difference of.")
                   } 
                   attrname = as.character(model.terms[[i]][[2]])
                       
                   S = c(S,
                         function(x,z){1*(x-z)})  
                   AlphaS_K = c(AlphaS_K, 1)
                   Snames = c(Snames, paste("WtoM_diff", attrname, sep="."))
               },
               MtoW_diff = {
                   
                   if(length(model.terms[[i]])<2){
                     stop("The 'MtoW_diff' term requires a first argument, the attribute to take the difference of.")
                   } 
                   attrname = as.character(model.terms[[i]][[2]])
                   if(length(model.terms[[i]])<3){
                     stop("The 'WtoM_diff' term requires a second argument, the difference level.")
                   } 
                   
                   S = c(S,
                         function(z,x){1*(z-x)})  
                   AlphaS_K = c(AlphaS_K, 1)
                   Snames = c(Snames, paste("MtoW_diff", attrname, sep="."))
               },
               # default
               {
                   message("Error: unrecognized rpm term: ", as.character(model.terms[[i]][2]))
               }
               
               )
    }

    names(S) = Snames
    names(Sr) = Snames
    names(Spaired) = Snames
    names(AlphaS_K) = Snames
    return(list(S=S,Sr=Sr,Spaired=Spaired,AlphaS_K=AlphaS_K))
}
