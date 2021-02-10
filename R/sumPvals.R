#' @title True Discovery Guarantee for p-Value Combinations
#' @description This function determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest, when using p-values as test statistics.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage sumPvals(G, S = seq(ncol(G)), alpha = 0.05, truncFrom = alpha, truncTo = max(alpha, 0.5),
#'          type = "vovk.wang", r = 0, nMax = 50)
#' @param G numeric matrix of p-values, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' The first transformation is the identity.
#' @param S vector of indices for the variables of interest.
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, p-values are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, p-values are not truncated.
#' @param type p-value combination (\code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{vovk.wang})
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param nMax maximum number of iterations.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Edgington: \code{-p}
#' \item Fisher: \code{-log(p)}
#' \item Pearson: \code{log(1-p)}
#' \item Liptak: \code{-qnorm(p)}
#' \item Cauchy: \code{tan(0.5 - p)/p}
#' \item Vovk and Wang: \code{- sign(r)p^r}
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1), where
#' \code{B} is the number of data transformations (rows in \code{G}).
#' @return \code{sumPvals} returns an object of class \code{sumSome}, containing
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
#' # generate matrix of p-values for 5 variables and 10 permutations
#' G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, seed = 42)
#' 
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#'  
#' # create object of class sumSome
#' # combination: harmonic mean (Vovk and Wang with r = -1)
#' res <- sumPvals(G, S, alpha = 0.4, r = -1)
#' 
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
#' @export

  
sumPvals <- function(G, S=seq(ncol(G)), alpha=0.05, truncFrom=alpha, truncTo=max(alpha, 0.5), type="vovk.wang", r=0, nMax=50){
  
  type = match.arg(tolower(type), c("fisher", "pearson", "liptak", "edgington", "cauchy", "vovk.wang"))
  res <- transf(G, truncFrom, truncTo, type, r)
  rm(G)
  
  out <- sumTest(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)
  return(out)
}