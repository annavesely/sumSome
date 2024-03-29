#' @title Simulating Matrix of Statistics
#' @description This function simulates a matrix of permutation statistics,
#' by performing a t test on normal data.
#' @usage simData(prop, m, B = 200, rho = 0, n = 50, alpha = 0.05, pw = 0.8, p = TRUE, seed = NULL)
#' @param prop proportion of non-null hypotheses.
#' @param m total number of variables.
#' @param B number of permutations, including the identity.
#' @param rho level of equicorrelation between pairs of variables.
#' @param n number of observations.
#' @param alpha significance level.
#' @param pw power of the t test.
#' @param p logical, \code{TRUE} to compute p-values, \code{FALSE} to compute t-scores.
#' @param seed seed.
#' @details The function applies the one-sample two-sided t test to a matrix of simulated data,
#' for \code{B} data permutations.
#' Data is obtained by simulating \code{n} independent observations from a multivariate normal distribution,
#' where a proportion \code{prop} of the variables has non-null mean.
#' This mean is such that the one-sample t test with significance level \code{alpha} has power equal to \code{pw}.
#' Each pair of distinct variables has equicorrelation \code{rho}.
#' @return \code{simData} returns a matrix where the \code{B} rows correspond to permutations (the first is the identity),
#' and the \code{m} columns correspond to variables.
#' The matrix contains p-values if \code{p} is \code{TRUE}, and t-scores otherwise.
#' The first columns (a proportion \code{prop}) correspond to non-null hypotheses.
#' @author Anna Vesely.
#' @examples
#' # generate matrix of p-values for 5 variables and 10 permutations
#' G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, seed = 42)
#' 
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#'  
#' # create object of class sumObj
#' # combination: harmonic mean (Vovk and Wang with r = -1)
#' res <- sumPvals(G, S, alpha = 0.4, r = -1)
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
#' @seealso
#' True discovery guarantee: \code{\link{sumStats}}, \code{\link{sumPvals}}
#' @importFrom pARI signTest
#' @export


simData <- function(prop, m, B=200, rho=0, n=50, alpha=0.05, pw=0.8, p=TRUE, seed=NULL){
  
  if(!is.numeric(prop) || !is.finite(prop)){stop("prop must be a number in [0,1]")}
  if(prop < 0 || prop > 1){stop("prop must be a number in [0,1]")}
  
  # check on m, B, n
  if(!is.numeric(m) || !is.finite(m) || m <= 0){stop("m must be a positive integer")}
  if(!is.numeric(B) || !is.finite(B) || B <= 0){stop("B must be a positive integer")}
  if(!is.numeric(n) || !is.finite(n) || n <= 0){stop("n must be a positive integer")}
  m <- round(m)
  B <- round(B)
  n <- round(n)
  
  # check on rho, alpha, pw
  if(!is.numeric(rho) || !is.finite(rho) || rho < -1 || rho > 1){stop("rho must be a number in [-1,1]")}
  if(!is.numeric(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  if(!is.numeric(pw) || !is.finite(pw) || pw <= 0 || pw >= 1){stop("pw must be a number in [0,1]")}
  
  if(!is.null(seed)){if(!is.numeric(seed) || !is.finite(seed)){stop("seed must be a finite integer")}}
  else{seed <- sample(seq(.Machine$integer.max), 1)}
  set.seed(round(seed))
  
  pi0 <- 1 - prop
  X <- pARI::simulateData(pi0=pi0, m=m, n=n, rho=rho, power=pw, alpha=alpha)
  
  res <- pARI::signTest(X, B=B, alternative="two.sided", seed=seed, rand=FALSE)
  
  if(p){G <- rbind(res$pv, t(res$pv_H0))}
  else{G <- abs(rbind(res$Test, t(res$Test_H0)))}
  return(G)
}



