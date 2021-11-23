#' @title True Discovery Guarantee for Generic Statistics
#' @description This function determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage sumStats(G, S = NULL, alternative = "greater", alpha = 0.05,
#'          truncFrom = NULL, truncTo = NULL, nMax = 50)
#' @param G numeric matrix of statistics, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' The first transformation is the identity.
#' @param S vector of indices for the variables of interest (if not specified, all variables).
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param nMax maximum number of iterations.
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1), where
#' \code{B} is the number of data transformations (rows in \code{G}).
#' @return \code{sumStats} returns an object of class \code{sumObj}, containing
#' \itemize{
#' \item \code{total}: total number of variables (columns in \code{G})
#' \item \code{size}: size of \code{S}
#' \item \code{alpha}: significance level
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{iterations}: number of iterations of the algorithm
#' }
#' @author Anna Vesely.
#' @examples
#' # generate matrix of t-scores for 5 variables and 10 permutations
#' G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, p = FALSE, seed = 42)
#'  
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#'  
#' # create object of class sumObj
#' res <- sumStats(G, S, alpha = 0.4, truncFrom = 0.7, truncTo = 0)
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
#' Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' True discovery guarantee using p-values: \code{\link{sumPvals}}
#' 
#' Access a \code{sumObj} object: \code{\link{discoveries}}, \code{\link{tdp}}, \code{\link{fdp}}
#' @export


sumStats <- function(G, S=NULL, alternative="greater", alpha=0.05, truncFrom=NULL, truncTo=NULL, nMax=50){
  
  if(is.null(S)){S <- seq(ncol(G))}
  
  alternative <- match.arg(tolower(alternative), c("greater", "lower", "two.sided"))
  res <- transf(G, truncFrom, truncTo, alternative, 1)
  rm(G)
  
  out <- sumTest(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)
  return(out)
}