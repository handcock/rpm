.recurse_summation <- function(x, sign){
  if(length(x)==1) {out <- list(x); attr(out,"sign")<-sign; out}
  else if(length(x)==2 && x[[1L]]=="+") .recurse_summation(x[[2L]],sign)
  else if(length(x)==2 && x[[1L]]=="-") .recurse_summation(x[[2L]],-sign)
  else if(length(x)==3 && x[[1L]]=="+") {
    l1 <- .recurse_summation(x[[2L]],sign)
    l2 <- .recurse_summation(x[[3L]],sign)
    out <- c(l1, l2)
    attr(out,"sign") <- c(attr(l1,"sign"), attr(l2,"sign"))
    out
  }
  else if(length(x)==3 && x[[1L]]=="-"){
    l1 <- .recurse_summation(x[[2L]],sign)
    l2 <- .recurse_summation(x[[3L]],-sign)
    out <- c(l1, l2)
    attr(out,"sign") <- c(attr(l1,"sign"), attr(l2,"sign"))
    out
  }
  else if(x[[1L]]=="(") .recurse_summation(x[[2L]], sign)
  else {out <- list(x); attr(out,"sign")<-sign; out}
}

#' Returns a list containing the terms in a given formula
#' @param object formula A formula having a right-hand-side that can be interpretated as a rpm specification.
#' returns a list containing terms in a given
#' formula, handling \code{+} and \code{-} operators and parentheses, and
#' keeping track of whether a term has a plus or a minus sign.
#'
#' @return
#' \code{list_rhs.formula} returns a list of formula terms, with an additional numerical vector attribute \code{"sign"} with of the same length, giving the corresponding term's sign as \code{+1} or \code{-1}.
#' @export
list_rhs.formula<-function(object){
  if (!is(object, "formula"))
    stop("Invalid formula of class ",sQuote(class(object)),".")
  
  .recurse_summation(ult(object), sign=+1)
}
#' @rdname ult
#'
#' @param value Replacement value for the `i`th element from the end.
#'
#' @note Due to the way in which assigning to a function is
#'   implemented in R, `ult(x) <- e` may be less efficient than
#'   `x[[length(x)]] <- e`.
#' 
#' @examples
#' (x <- c(1:5))
#' (ult(x) <- 6)
#' (ult(x, 2) <- 7) # 2nd last.
#' x
#'
#' \dontshow{
#' stopifnot(all(x == c(1:3, 7L, 6L)))
#' }
#'
#' @export
`ult<-` <- function(x, i=1L, value){
  x[[length(x)-i+1L]] <- value
  x
}
#' Extract or replace the *ult*imate (last) element of a vector or a list, or an element counting from the end.
#'
#' @param x a vector or a list.
#' @param i index from the end of the list to extract or replace (where 1 is the last element, 2 is the penultimate element, etc.).
#'
#' @return An element of `x`.
#'
#' @examples
#' x <- 1:5
#' (last <- ult(x))
#' (penultimate <- ult(x, 2)) # 2nd last.
#'
#' \dontshow{
#' stopifnot(last==5)
#' stopifnot(penultimate==4)
#' }
#'
#' @export
ult <- function(x, i=1L){
  x[[length(x)-i+1L]]
}
