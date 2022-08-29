#' Terms used in a Revealed Preference Matchings Model
#' 
#' The function \code{\link{rpm}} is used to fit a revealed preference model
#' for men and women of certain
#' characteristics (or shared characteristics) of people of the opposite sex.
#' The model assumes a one-to-one stable matching using an observed set of
#' matchings and a set of (possibly dyadic) covariates to 
#' estimate the parameters for
#' linear equations of utilities.
#' It does this using an approximate likelihood based on ideas from Menzel (2015).
#' 
#' The model represents the dyadic utility functions as deterministic linear utility functions of
#' dyadic variables. These utility functions are functions of observed characteristics of the women
#' and men.
#' These functions are entered as terms in the function call
#' to \code{\link{rpm}}.  This page describes the possible terms (and hence
#' linear utility functions) included in \code{\link[=rpm-package]{rpm}} package.
#'
#' @docType methods
#' @name rpm-terms
#' @aliases rpm-terms rpm.terms terms-rpm terms.rpm absdiff
#' W_atleast W_atmost W_cov W_factor W_greaterthan 
#' M_atleast M_atmost M_cov M_factor M_greaterthan 
#' WtoM_diff MtoW_diff
#' match 
#' diff
#' mix 
#' 
#' @section Specifying models:
#' 
#' Terms to \code{\link{rpm}} are specified by a formula to represent the
#' pairings and covariates This is done via a \code{formula}, that is,
#' an formula object, of the form \code{~ <term 1> + <term 2> ...}, where
#' \code{<term 1>}, \code{<term 2>}, etc, are each terms chosen
#' from the list given below.
#' \describe{
#'  \item{\code{absdiff(attr)} (quantitative attribute),
#'    \code{absdiff(attr)} (quantitative attribute)}{\emph{Absolute difference:}
#'    The \code{attr} argument specifies a quantitative attribute 
#'    This term adds one statistic to the model equaling
#'    \code{abs(attr[i]-attr[j])} for all women-man dyad (i,j).
#'   }
#'  \item{\code{W_greaterthan(attr)}}{\emph{Women's value
#'    greater than the men's value} Adds one statistic
#'    indicating if the women's value exceeds the men's value.
#'   }
#'  \item{\code{M_greaterthan(attr)}}{\emph{Men's value
#'    greater than the women's value} Adds one statistic
#'    indicating if the men's value exceeds the women's value.
#'   }
#'  \item{\code{W_atleast(attr,threshold=0)}}{\emph{Values
#'    greater than or equal to a threshold} Adds one statistic
#'    indicating if the women's value of the attribute equals or exceeds
#'    \code{threshold}.
#'   }
#'  \item{\code{W_atmost(threshold=0)}}{\emph{Values
#'    less than or equal to a threshold} Adds one statistic
#'    indicating if the women's value equals or is exceeded by
#'    \code{threshold}.
#'   }
#'  \item{\code{W_cov(attr)} (quantitative attribute),
#'    \code{W_cov(attr)} (quantitative attribute)
#'    }{\emph{Main effect of a covariate for women:}
#'    The \code{attr} argument specifies a quantitative attribute
#'    This term adds a single statistic equaling the 
#'    value of \code{attr(i)} for women \eqn{i} in the dyad.
#'    For categorical attributes,
#'    see \code{W_factor}.
#'   }
#'  \item{\code{diff(attr)} (quantitative attribute), \code{diff(attr)}
#'    (quantitative attribute)}{\emph{Woman's Gap:}
#'    The \code{attr} argument specifies a quantitative attribute
#'    This term adds one statistic to the model
#'    being \code{attr[i]-attr[j]} for women \eqn{i} and man \eqn{j}.
#'    Specifically, it is the excess of the woman's value over the man's value.
#'   }
#'  \item{\code{WtoM_diff(attr, diff)} (ordinal categorical attribute), \code{WtoM_diff(attr)}
#'    (ordinal categorical discrete attribute)}{\emph{Woman's Gap:}
#'    The \code{attr} argument specifies a ordinal categorical attribute
#'    This term adds one statistic to the model
#'    being an indicator that \code{attr[i]=attr[j]+diff} for women \eqn{i} and man \eqn{j}.
#'    Specifically, it indicates if the woman's value is \code{diff} higher than the man's value.
#'   }
#'  \item{\code{MtoW_diff(attr, diff)} (ordinal categorical attribute), \code{MtoW_diff(attr)}
#'    (ordinal categorical discrete attribute)}{\emph{Man's Gap:}
#'    The \code{attr} argument specifies a ordinal categorical attribute
#'    This term adds one statistic to the model
#'    being an indicator that \code{attr[j]=attr[i]+diff} for women \eqn{i} and man \eqn{j}.
#'    Specifically, it indicates if the man's value is \code{diff} higher than the woman's value.
#'   }
#'  \item{\code{MtoW_diff(attr)} (quantitative attribute), \code{MtoW_diff(attr)}
#'    (quantitative attribute)}{\emph{Difference:}
#'    The \code{attr} argument specifies a quantitative attribute
#'    This term adds one statistic to the model
#'    \code{attr[j]-attr[i]} for women \eqn{i} and man \eqn{j}.
#'   }
#'  \item{\code{W_factor(attr, base=1, levels=-1)} (categorical attribute),
#'    \code{W_factor(attr, base=1, levels=-1)} (categorical attribute)
#'    }{\emph{Factor attribute effect for women:}
#'    The \code{attr} argument specifies a categorical attribute
#'    This term adds
#'    multiple statistics to the model, one for each of (a subset of) the
#'    unique values of the \code{attr} attribute. Each of these statistics
#'    indicates if the women's has that attribute.
#'   }
#'  \item{\code{match(attr, diff=FALSE, collapse=NULL)}
#'    }{\emph{Attribute-based homophily effect:}
#'    The \code{attr} argument specifies a categorical attribute
#'    This term adds one statistic to the model
#'    unless \code{diff} is set to \code{TRUE}, in which case the term adds multiple 
#'    statistics to the model, one for each of (a subset of) the unique values of the \code{attr}
#'    attribute. 
#'    If \code{diff} is set to \code{TRUE}, the optional argument \code{collapse} control what dyads 
#' are collapsed (or pooled).
#' Specifically, it is a list of indices of attribute values which are to be collapsed into a
#' single term. For example, \code{collapse=list(c(1,4))} will collapse the \code{(1,1)} and the 
#' \code{(4,4)} dyads into a single term (and group). Multiple lists can be included with arbitrary numbers of 
#'  dyads in a group. 
#'   }
#'  \item{\code{mix(attr, base=NULL, collapse=NULL)}
#'    }{\emph{Attribute mixing:} The \code{attr} argument specifies a categorical attributes
#'    By default, this term adds one statistic to
#'    the model for each possible pairing of attribute values. The
#'    statistic indicates if the dyad 
#'    has that pairing of values.
#'    In other words, this term produces one statistic for
#'    every entry in the mixing matrix for the attribute(s). The ordering of
#'    the attribute values is lexicographic: alphabetical (for nominal categories) or
#'    numerical (for ordered categories).
#'    The optional argument \code{base} control what statistics are 
#'    included in the model, specifically it lists the index of the omitted terms (in order).
#'    For example, \code{base=2} omits the second term. 
#'    The optional argument \code{collapse} control what dyads are collapsed (or pooled).
#' Specifically, it is a list of lists. Each element of the list is a list of dyads which are to be collapsed into a
#' single term. For example, \code{collapse=list(list(c(1,4),c(2,4)))} will collapse the \code{(1,4)} and the 
#' \code{(2,4)} dyads into a single term (and group). Multiple lists can be included with arbitrary numbers of 
#'  dyads in a group. 
#'   }
#' }
#' @seealso \code{\link[=rpm-package]{rpm}} package,
#' \code{\link{rpm}}
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' 
#' @keywords models
#' @examples
#' library(rpm)
#' data(fauxmatching)
#' fit <- rpm(~match("edu") + WtoM_diff("edu",3),
#'           Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
#'           X_w="X_w", Z_w="Z_w",
#'           pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
#'           sampled="sampled")
#' summary(fit)
NULL
