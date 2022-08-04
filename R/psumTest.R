#' @title True Discovery Guarantee - Parametric
#' @description Internal function.
#' It determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage psumTest(g, S, alpha, cvs)
#' @param g numeric vector of statistics.
#' @param S vector of indices for the variables of interest.
#' @param alpha significance level.
#' @param cvs numeric vector of critical values for summed statistics considering \code{1:m} hypotheses.
#' @return \code{sumTestPar} returns an object of class \code{sumObj}, containing
#' \itemize{
#' \item \code{total}: total number of variables (length of \code{g})
#' \item \code{size}: size of \code{S}
#' \item \code{alpha}: significance level
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{iterations}: number of iterations of the algorithm
#' }
#' @author Xu Chen, Anna Vesely.
#' @noRd


psumTest <- function(g, S, alpha, cvs){
  
  m <- length(g)
  s <- length(S)
  
  notS <- (1:m)[-S]
  u0 <- max(g)
  v0 <- min(g)
  u <- c(v0, sort(g[S]), rep(u0+1, m-s+1))
  v <- c(v0, sort(g[notS]), rep(u0+1, s+1))
  
  td <- findDiscSum(s, m, u, v, cvs)
  out <- sumObj(m, s, alpha, td, td, NULL)
  return(out)
}
