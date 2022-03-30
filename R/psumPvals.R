#' @title True Discovery Guarantee for p-Value Combinations - Parametric
#' @description This function determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest, when using p-values as test statistics.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage psumPvals(g, S = NULL, alpha = 0.05, type = "harmonic.dep")
#' @param g numeric vector of p-values.
#' @param S vector of indices for the variables of interest (if not specified, all variables).
#' @param alpha significance level.
#' @param type p-value combination among \code{harmonic.dep}, \code{harmonic.ind}, \code{fisher.ind}, \code{cauchy.ind}
#' (see details).
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Harmonic mean: \code{1/p}
#' \item Fisher: \code{-2log(p)}
#' \item Cauchy: \code{tan[(0.5 - p)pi]} with \code{pi}=3.142
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details The \code{type} determines the vector of critical values as following.
#' \itemize{
#' \item Harmonic mean under dependence: valid under general dependence (Vovk and Wang, 2020)
#' \item Harmonic mean under independence: valid under independence, anti-conservative under general dependence (Wilson, 2019)
#' \item Fisher under independence:  valid under independence, anti-conservative under general dependence (Fisher, 1925)
#' \item Cauchy: valid under independence and perfect dependence, approximately valid under general dependence (Liu and Xie, 2020)
#' }
#' @return \code{sumPvalsPar} returns an object of class \code{sumObj}, containing
#' \itemize{
#' \item \code{total}: total number of variables (length of \code{g})
#' \item \code{size}: size of \code{S}
#' \item \code{alpha}: significance level
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{iterations}: number of iterations of the algorithm
#' }
#' @author Xu Chen, Jelle Goeman.
#' @examples
#' # generate vector of p-values for 5 variables
#' g <- as.vector(simData(prop = 0.6, m = 5, B = 1, alpha = 0.4, seed = 42))
#' 
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#'  
#' # create object of class sumObj
#' # combination: harmonic mean under general dependence
#' res <- psumPvals(g, S, alpha = 0.4, type = "harmonic.dep")
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
#' Tian, J., Chen, X., Katsevich, E., Goeman, J. J., and Ramdas, A. (2021). Large-scale simultaneous inference under dependence. Pre-print arXiv:2102.11253.
#' 
#' @seealso
#' True discovery guarantee using generic statistics (parametric): \code{\link{psumStats}}
#' 
#' Access a \code{sumObj} object: \code{\link{discoveries}}, \code{\link{tdp}}, \code{\link{fdp}}
#' @export


psumPvals <- function(g, S=NULL, alpha=0.05, type="harmonic.dep"){
  
  if(!is.vector(g) || !is.numeric(g) || !all(is.finite(g))){stop("g must be a vector of finite numbers")}
  if(length(g)==0){stop("g must be a vector of finite numbers")}
  if(!all(g >= 0) || !all(g <= 1)){stop("g must be a vector of pvalues")}
  
  if(is.null(S)){S <- seq(length(g))}
  
  type = match.arg(tolower(type), c("harmonic.dep", "harmonic.ind", "fisher.ind", "cauchy"))
  
  if (type=="harmonic.dep" || type=="harmonic.ind"){g <- 1/g}
  else if(type=="fisher.ind"){g <- -2*log(g)}
  else if(type=="cauchy") {g <- tan((0.5 - g)*pi)}
  
  if(!all(is.finite(g))){stop("Transformation produced infinite values")}
  
  cvs <- generateCV(length(g), type, alpha)
  
  out <- psumTest(g, S, alpha, cvs)
  return(out)
}