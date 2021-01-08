#' @title True Discovery Guarantee for p-Value Combinations
#' @description This function determines a lower confidence bound for the number of true discoveries within a set of interest,
#' when using p-values as test statistics. The bound remains valid under post-hoc selection.
#' @usage sumSomePvals(G, S, alpha = 0.05, truncFrom = alpha, truncTo = 1, type = "vovk.wang", r = 0,
#'              nMax = 10000)
#' @param G numeric matrix of p-values, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' The first transformation is the identity.
#' @param S vector of indices for the variables of interest.
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, p-values are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, it is set to 0.5 for Pearson's and Liptak's methods, and 1 in the other cases. 
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
#' @return \code{sumSomePvals} returns a list containing \code{summary} (vector) and \code{iterations} (number of iterations).
#' The vector \code{summary} contains:
#' \itemize{
#' \item \code{size}: size of \code{S}
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the true discovery proportion in \code{S}
#' \item \code{maxTD}: maximum value of \code{TDP} that could be found under convergence of the algorithm
#' }
#' @author Anna Vesely.
#' @examples
#' G <- matrix(
#'  c(0.05, 0.20, 0.24, 0.45, 0.54,
#'    0.47, 0.34, 0.53, 0.99, 0.18,
#'    0.14, 0.21, 0.98, 0.32, 0.47,
#'    0.19, 0.45, 0.85, 0.50, 0.92,
#'    0.99, 0.09, 0.52, 0.39, 0.37,
#'    0.87, 0.89, 0.44, 0.38, 0.55),
#'  ncol=5, byrow=TRUE)
#'  
#' # harmonic mean (Vovk and Wang with r=-1)
#' sumSomePvals(G, S = c(1,2), alpha = 0.4, r = -1)
#' @export

  
sumPvals <- function(G, S, alpha=0.05, truncFrom=alpha, truncTo=0.5, type="vovk.wang", r=0, nMax=10000){
  
  type = match.arg(tolower(type), c("fisher", "pearson", "liptak", "edgington", "cauchy", "vovk.wang"))
  res <- transf(G, truncFrom, truncTo, type, r)
  rm(G)
  
  out <- sumTest(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)
  return(out)
}