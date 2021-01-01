#' @title True Discovery Guarantee for p-Values
#' @description This function determines a lower confidence bound for the number of true discoveries within a set of interest. The bound remains valid under post-hoc selection.
#' @usage sumSome.pvalues(G, S, alpha = 0.05, truncFrom = alpha, truncTo = 1, type = "fisher", r = 1, nMax = 10000)
#' @param G numeric matrix of p-values, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' @param S vector of indices for the variables of interest.
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, p-values are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, p-values are not truncated.
#' @param nMax maximum number of iterations.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' @details Pearson's and Liptak's transformations produce infinite values in \code{1}.
#' For such transformations, \code{truncTo} is coerced to be not greater than \code{1 -  .Machine$double.eps}.
#' @details An error message is returned if the p-value transformation produces infinite values.
#' @return \code{sumSome.pvalues} returns a list containing \code{summary} (vector) and \code{iterations} (number of iterations).
#' The vector \code{summary} contains
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
#'  c(0.05, 0.12, 0.21, 0.66, 0.66,
#'   0.66, 0.50, 0.66, 0.79, 0.21,
#'   0.01, 0.34, 0.79, 0.50, 0.66,
#'   0.01, 0.66, 0.79, 0.66, 0.79,
#'   0.79, 0.05, 0.66, 0.66, 0.50,
#'   0.02, 0.79, 0.66, 0.50, 0.66),
#' ncol=5, byrow=TRUE)
#' S <- c(1,2)
#' sumSome.pvalues(G, S, alpha = 0.4, type = "fisher")
#' @export


sumSome.pvalues <- function(G, S, alpha=0.05, truncFrom=alpha, truncTo=1,
                            type="fisher", r=1, nMax=10000){
  
  type = match.arg(tolower(type), c("fisher", "pearson", "liptak", "edgington", "cauchy", "vovk.wang"))
  res <- transf(G, truncFrom, truncTo, type, r)
  rm(G)
  out <- sumSome.internal(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)
  return(out)
}