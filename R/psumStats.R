#' @title True Discovery Guarantee for Generic Statistics - Parametric
#' @description This function determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage psumStats(g, S = NULL, alpha = 0.05, cvs)
#' @param g numeric vector of statistics.
#' @param S vector of indices for the variables of interest (if not specified, all variables).
#' @param alpha significance level.
#' @param cvs numeric vector of critical values for summed statistics considering \code{1:m} hypotheses.
#' @return \code{sumStatsPar} returns an object of class \code{sumObj}, containing
#' \itemize{
#' \item \code{total}: total number of variables (length of \code{g})
#' \item \code{size}: size of \code{S}
#' \item \code{alpha}: significance level
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{iterations}: number of iterations of the algorithm
#' }
#' @author Xu Chen, Anna Vesely.
#' @examples
#' # generate vector of statistics for 5 variables (Fisher transformation of p-values)
#' g <- as.vector(simData(prop = 0.6, m = 5, B = 1, alpha = 0.4, seed = 42))
#' g <- -2 * log(g)
#' 
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#' 
#' # vector of critical values
#' cvs <- qchisq(p = 0.4, df = 2 * seq(5), lower.tail=FALSE)
#'  
#' # create object of class sumObj
#' res <- psumStats(g, S, alpha = 0.4, cvs = cvs)
#' res
#' summary(res)
#' 
#' # lower confidence bound for the number of true discoveries in S
#' discoveries(res)
#' 
#' # lower confidence bound for the true discovery proportion in S
#' tdp(res)
#' 
#' # upper confidence bound for the false discovery proportion in S
#' fdp(res)
#' @references
#' Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.
#' 
#' Tian, J., Chen, X., Katsevich, E., Goeman, J. J., and Ramdas, A. (2021). Large-scale simultaneous inference under dependence. Scandinavian Journal of Statistics, to appear. (Pre-print arXiv:2102.11253)
#' 
#' @seealso
#' True discovery guarantee using p-values (parametric): \code{\link{psumPvals}}
#' 
#' Access a \code{sumObj} object: \code{\link{discoveries}}, \code{\link{tdp}}, \code{\link{fdp}}
#' @export


psumStats <- function(g, S=NULL, alpha=0.05, cvs){
  
  if(!is.vector(g) || !is.numeric(g) || !all(is.finite(g))){stop("g must be a vector of finite numbers")}
  if(length(g)==0){stop("g must be a vector of finite numbers")}
  
  if(is.null(S)){S <- seq(length(g))}
  if(!is.vector(S) || !is.numeric(S) || !all(is.finite(S))){stop("S must be a vector of finite integers")}
  if(!all(floor(S)==S)){stop("S must be a vector of finite integers")}
  if(!all(S >= 0) || !all(S <= length(g))){stop("S must contain indices between 1 and the total number of variables")}
  S <- unique(S)
  
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  
  if(!is.vector(cvs) || !is.numeric(cvs) || !all(is.finite(cvs))){stop("cvs must be a vector of finite numbers")}
  if(length(g) != length(cvs)){stop("Invalid dimensions: g and cvs must have the same length")}
  
  out <- psumTest(g, S, alpha, cvs)
  return(out)
}