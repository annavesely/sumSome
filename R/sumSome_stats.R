#' @title True Discovery Guarantee
#' @description This function determines a lower confidence bound for the number of true discoveries within a set of interest.
#' The bound remains valid under post-hoc selection.
#' @usage sumSome.stats(G, S, alternative = "greater", alpha = 0.05, truncFrom = NULL, truncTo = NULL, nMax = 10000)
#' @param G numeric matrix of statistics, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' @param S vector of indices for the variables of interest.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param nMax maximum number of iterations.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' @return \code{sumSome.stats} returns a list containing \code{summary} (vector) and \code{iterations} (number of iterations).
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
#'  c(6,5,4,1,1,
#'   1,2,1,0,4,
#'   8,3,0,2,1,
#'   8,1,0,1,0,
#'   0,6,1,1,2,
#'   7,0,1,2,1),
#'  ncol=5, byrow=TRUE)
#' S <- c(1,2)
#' sumSome.stats(G, S, alternative = "greater", alpha = 0.4, truncFrom = 2, truncTo = 0)
#' @export


sumSome.stats <- function(G, S, alternative="greater", alpha=0.05, truncFrom=NULL, truncTo=NULL, nMax=10000){
  
  alternative <- match.arg(tolower(alternative), c("greater", "lower", "two.sided"))
  res <- transf(G, truncFrom, truncTo, alternative, 1)
  rm(G)
  out <- sumSome.internal(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)
  return(out)
}