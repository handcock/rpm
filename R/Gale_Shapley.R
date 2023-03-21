#' This is the version of Gale-Shapley stable matching algorithm (translated from the Matlab code in Menzel (2015)).
#' 
#' This code allows the self-matched option
#' 
#' @param U The utility matrix for the women's side. Each row is a woman, each column is a man.
#' The matrix entry (i,j) is the utility that woman \code{i} gains from pairing with man \code{j}. 
#' In other words, the utility is computed from woman \code{i}'s perspective.
#' @param V The utility matrix for the men's side. Each column is a man, each row is a woman.
#' The matrix entry (i,j) is the utility that man \code{j} gains from pairing with woman \code{i}. 
#' In other words, the utility is computed from man \code{j}'s perspective.
#' @param return.data.frame logical Should a \code{data.frame} of the matching be returned instead of the
#' paring matrix mu?
#' @param cpp logical Should the \code{Rcpp} version of the code be used. This is much faster and uses a lot less memory.
#' @param nmax count The maximum number of iterations of the inner loop within the Gale-Shapley algorithm.
#' This can be reduced to speed up the algorithm at the potential cost of many partnerships being non-equilibruim.
#' @return The function return depends on the \code{return.data.frame} value. 
#' If TRUE, it returns
#' \item{data.frame}{a two-column \code{data.frame} with the first column a women's index and the second column the 
#' men's index of their partner. It has as many rows as there are partnerships.}
#' If FALSE, it returns the following matrix: 
#' \item{mu}{If \code{cpp=TRUE}, a vector of length the number of women (\code{nrow(U)}) with the 
#' index of the matching man (i.e., the index is the row in \code{V} of the man). If there is no
#' matching man, the index is 0. This can be used to reconstruct the matching matrix.
#' If \code{cpp=FALSE}, the matching matrix, where 1 represents a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows is the same as the 
#' rows in \code{U}. The order of the columns is the same as the columns in \code{V}.}
#' @seealso rpm
#' @references Goyal, Handcock, Jackson. Rendall and Yeung (2023).
#' \emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation}
#' \emph{Journal of the Royal Statistical Society}, A. \doi{10.18637/jss.v024.i07} 
#' Menzel, K. (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @export Gale_Shapley
#'@importFrom Rcpp evalCpp
#'@useDynLib rpm

Gale_Shapley <- function(U,V,return.data.frame=FALSE, cpp=TRUE, nmax=10*nrow(U)){

# women (U) are the proposing side, with women listed as rows

 nw <- nrow(U)
 nm <- ncol(U)-1

 if(cpp){
    # These are the integer versions
    Ua <- matrixStats::rowRanks(U)
    U <- sweep(Ua[,-1],1,Ua[,1],"-")
    O=t(apply(-U,1,order))
    Ua <- matrixStats::rowRanks(-U[,-1])
    V <- matrixStats::rowRanks(V)
    V <- sweep(V[,-1],1,V[,1],"-")
    storage.mode(U) <- "integer"
    storage.mode(V) <- "integer"
    storage.mode(O) <- "integer"
    storage.mode(Ua) <- "integer"
    mu <- GSi_NTU_O(U,V,O,Ua,nmax=nmax)
    rm(U,V,O,Ua)
 }else{
   U <- sweep(U[,-1],1,U[,1],"-")
   V <- sweep(V[,-1],1,V[,1],"-")
   
   for (i in 1:nmax){
     
     # check dimensions!
     # no need to worry about ties
     Prop <- PropMax(U)
     PropV <- Prop*t(V)
     Rej <- (PropV < sweep(Prop,2,colMax(rbind(PropV,0)),"*"))
     U[Rej&Prop] = -1
     
     if (sum(sum(Rej))==0){
       break
     }
   }
   
   mu <- sweep(U,1,rowMax(cbind(U,0)),"==")
 }
 
 if(return.data.frame){
   if(cpp){
     w_cpp <- as.vector(mu)
     m_cpp <- rep(0,nm)
     m_cpp[w_cpp[w_cpp>0]] <- (1:nw)[w_cpp>0] 
     w_cpp[w_cpp>0] <- w_cpp[w_cpp>0] + nw
     pair_id_int <- c(w_cpp, m_cpp)
     pair_id=as.character(as.integer(pair_id_int))
     pair_id[pair_id=="0"] <- NA
   }else{
     pair_id_w=rep(NA,length=nw)
     pair_id_w[unlist(apply(mu,1,function(x){any(x>0.5)}))] <- (nw+(1:nm))[unlist(apply(mu,1,function(x){which(x>0.5)}))]
     pair_id_m=rep(NA,length=nm)
     pair_id_m[unlist(apply(mu,2,function(x){any(x>0.5)}))] <- (1:nw)[unlist(apply(mu,2,function(x){which(x>0.5)}))]
     
     pair_id_int=c(pair_id_w,pair_id_m)
     pair_id=as.character(as.integer(pair_id_int))
     pair_id[is.na(pair_id_int)] <- NA
   }
   mu <- data.frame(pid=as.character(as.integer(1:(nw+nm))),
                    gender=rep(c("F","M"),c(nw,nm)),
                    pair_id=pair_id)
   return(mu)
 }else{
   return(1*mu) # for matrix return
 }
}
