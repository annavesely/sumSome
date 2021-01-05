#' @title Simulating Matrix of Statistics
#' @description This function simulates a matrix of statistics that may be used
#' to test \code{\link{sumSome}} and \code{\link{sumSome.pvals}}.
#' @usage simulateMat(prop, m = 1000, B = 200, rho = 0, n = 50, alpha = 0.05, power = 0.8, pvalues = TRUE,
#'             rand = FALSE, seed = NULL)
#' @param prop proportion of non-null hypotheses.
#' @param m total number of variables.
#' @param B number of permutations, including the identity.
#' @param rho level of equicorrelation between pairs of variables.
#' @param n number of observations.
#' @param alpha significance level.
#' @param power power of the one-sample t test.
#' @param pvalues logical, \code{TRUE} to compute p-values, \code{FALSE} to return t-scores.
#' @param rand logical, \code{TRUE} to compute p-values by permutation distribution.
#' @param seed seed.
#' @details The function applies the one-sample two-sided t test to a matrix of simulated data,
#' for \code{B} data permutations.
#' Data is obtained by simulating \code{n} independent observations from a multivariate normal distribution,
#' where a proportion \code{prop} of the variables has non-null mean.
#' This mean is such that the one-sample t test with significance level \code{alpha} has power equal to \code{power}.
#' Each pair of distinct variables has equicorrelation \code{rho}.
#' @return \code{simulateMat} returns a matrix where the \code{B} rows correspond to permutations (the first is the identity),
#' and the \code{m} columns correspond to variables.
#' The matrix contains p-values if \code{pvalues} is \code{TRUE}, and t-scores otherwise.
#' The first columns (a proportion \code{prop}) correspond to non-null hypotheses.
#' @author Anna Vesely.
#' @examples
#' # simulate 10 x 5 matrix of p-values
#' # the first 3 variables correspond to non-null hypotheses
#' G <- simulateMat(prop = 0.6, m = 5, B = 10, alpha = 0.4, seed = 42)
#' 
#' sumSome.pvals(G, S = c(1,2,3), alpha = 0.4, type = "vovk.wang", r = -1)
#' @export


simulateMat <- function(prop, m=1000, B=200, rho=0, n=50, alpha=0.05, power=0.8, pvalues=TRUE, rand=FALSE, seed=NULL){
  
  if(!is.numeric(prop) || !is.finite(prop)){stop("prop must be a number in (0,1)")}
  if(prop < 0 || prop > 1){stop("prop must be a number in [0,1]")}
  
  # check on m, B, n
  if(!is.numeric(m) || !is.finite(m) || m <= 0){stop("m must be a positive integer")}
  if(!is.numeric(B) || !is.finite(B) || B <= 0){stop("B must be a positive integer")}
  if(!is.numeric(n) || !is.finite(m) || n <= 0){stop("n must be a finite integer")}
  m <- round(m)
  B <- round(B)
  n <- round(n)
  
  # check on rho, alpha, power
  if(!is.numeric(rho) || !is.finite(rho)){stop("rho must be a number in [-1,1]")}
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(!is.numeric(power) || !is.finite(power)){stop("power must be a number in [0,1]")}
  if(rho < -1 || rho > 1){stop("rho must be a number in [-1,1]")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  if(power < 0 || power > 1){stop("power must be a number in [0,1]")}
  
  if(!is.null(seed)){
    if(!is.numeric(seed) || !is.finite(seed)){stop("seed must be a finite integer")}
    set.seed(round(seed))
  }
  
  pi0 <- 1 - prop
  X <- pARI::simulateData(pi0=pi0, m=m, n=n, rho=rho, power=power, alpha=alpha)
  
  res <- pARI::signTest(X, B=B, alternative="two.sided", seed=seed, rand=rand)
  
  if(pvalues){G <- rbind(res$pv, t(res$pv_H0))}
  else{G <- rbind(res$Test, t(res$Test_H0))}
  return(G)
}



