% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gale_Shapley.R
\name{Gale_Shapley}
\alias{Gale_Shapley}
\title{This is the version of Gale-Shapley stable matching algorithm (translated from the Matlab code in Menzel (2015)).}
\usage{
Gale_Shapley(U, V, return.data.frame = FALSE, cpp = TRUE, nmax = 10 * nrow(U))
}
\arguments{
\item{U}{The utility matrix for the women's side. Each row is a woman, each column is a man.
The matrix entry (i,j) is the utility that woman \code{i} gains from pairing with man \code{j}. 
In other words, the utility is computed from woman \code{i}'s perspective.}

\item{V}{The utility matrix for the men's side. Each column is a man, each row is a woman.
The matrix entry (i,j) is the utility that man \code{j} gains from pairing with woman \code{i}. 
In other words, the utility is computed from man \code{j}'s perspective.}

\item{return.data.frame}{logical Should a \code{data.frame} of the matching be returned instead of the
paring matrix mu?}

\item{cpp}{logical Should the \code{Rcpp} version of the code be used. This is much faster and uses a lot less memory.}

\item{nmax}{count The maximum number of iterations of the inner loop within the Gale-Shapley algorithm.
This can be reduced to speed up the algorithm at the potential cost of many partnerships being non-equilibruim.}
}
\value{
The function return depends on the \code{return.data.frame} value. 
If TRUE, it returns
\item{data.frame}{a two-column \code{data.frame} with the first column a women's index and the second column the 
men's index of their partner. It has as many rows as there are partnerships.}
If FALSE, it returns the following matrix: 
\item{mu}{If \code{cpp=TRUE}, a vector of length the number of women (\code{nrow(U)}) with the 
index of the matching man (i.e., the index is the row in \code{V} of the man). If there is no
matching man, the index is 0. This can be used to reconstruct the matching matrix.
If \code{cpp=FALSE}, the matching matrix, where 1 represents a pairing, 0 otherwise. 
Each row is a woman, each column is a man. The order of the rows is the same as the 
rows in \code{U}. The order of the columns is the same as the columns in \code{V}.}
}
\description{
This code allows the self-matched option
}
\references{
Goyal, Shuchi; Handcock, Mark S.; Jackson, Heide M.; Rendall, Michael S. and Yeung, Fiona C. (2023).
\emph{A Practical Revealed Preference Model for Separating Preferences and Availability Effects in Marriage Formation},
\emph{Journal of the Royal Statistical Society}, A. \doi{10.1093/jrsssa/qnad031} 

Dagsvik, John K. (2000) \emph{Aggregation in Matching Markets} \emph{International Economic Review}, Vol. 41, 27-57.
JSTOR: https://www.jstor.org/stable/2648822, \doi{10.1111/1468-2354.00054}

Menzel, Konrad (2015).
\emph{Large Matching Markets as Two-Sided Demand Systems}
Econometrica, Vol. 83, No. 3 (May, 2015), 897-941. \doi{10.3982/ECTA12299}
}
\seealso{
rpm
}
\keyword{models}
