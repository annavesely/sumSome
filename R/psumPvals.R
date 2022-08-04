#' @title True Discovery Guarantee for p-Value Combinations - Parametric
#' @description This function determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest, when using p-values as test statistics.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage psumPvals(g, S = NULL, alpha = 0.05, type = "harmonic.dep", r = NULL)
#' @param g numeric vector of p-values.
#' @param S vector of indices for the variables of interest (if not specified, all variables).
#' @param alpha significance level.
#' @param type p-value combination among \code{fisher}, \code{fisher.dep}, \code{pearson}, \code{liptak}, \code{cauchy}, \code{vovk.wang}, \code{harmonic.dep}, \code{harmonic.ind}
#' (see details).
#' @param r parameter for Vovk and Wang's p-value combination for general dependence.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Fisher: \code{-2log(p)}
#' \item Pearson: \code{2log(1-p)}
#' \item Liptak: \code{qnorm(1-p)}
#' \item Cauchy: \code{tan[(0.5-p)pi]} with \code{pi=3.142}
#' \item P-value combination for general dependence: \code{-sign(r)*p^r} for \code{r!=0} & \code{-log(p)} for \code{r=0}
#' \item Harmonic mean: \code{1/p}
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details The \code{type} determines the vector of critical values as following.
#' \itemize{
#' \item Fisher (independence, \code{fisher}): valid under independence, anti-conservative otherwise (Fisher, 1925)
#' \item Fisher (dependence, \code{fisher.dep}): valid under general dependence (Vovk and Wang, 2020)
#' \item Pearson (\code{pearson}): valid under independence, anti-conservative otherwise (Pearson, 1933)
#' \item Liptak (\code{liptak}): valid under independence, conservative or anti-conservative otherwise depending on the dependence structure among the tests (Liptak, 1958; unweighted version, same as Stouffer’s method (Stouffer et al., 1949))
#' \item Cauchy (\code{cauchy}): valid under independence and perfect dependence, approximately valid otherwise under bivariate normality assumption (Liu and Xie, 2020)
#' \item P-value combination for general dependence (\code{vovk.wang}): valid under general dependence (Vovk and Wang, 2020)
#' \item Harmonic mean (dependence, \code{harmonic.dep}): valid under general dependence (Vovk and Wang, 2020)
#' \item Harmonic mean (independence, \code{harmonic.ind}): valid under independence, anti-conservative otherwise (Wilson, 2019)
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
#' @author Xu Chen, Anna Vesely.
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
#' Tian, J., Chen, X., Katsevich, E., Goeman, J. J., and Ramdas, A. (2021). Large-scale simultaneous inference under dependence. Scandinavian Journal of Statistics, to appear. (Pre-print arXiv:2102.11253)
#' 
#' @seealso
#' True discovery guarantee using generic statistics (parametric): \code{\link{psumStats}}
#' 
#' Access a \code{sumObj} object: \code{\link{discoveries}}, \code{\link{tdp}}, \code{\link{fdp}}
#' @export

psumPvals <- function(g, S=NULL, alpha=0.05, type="harmonic.dep", r=NULL){
  
  if(!is.vector(g) || !is.numeric(g) || !all(is.finite(g))){stop("g must be a vector of finite numbers")}
  if(length(g)==0){stop("g must be a vector of finite numbers")}
  if(!all(g >= 0) || !all(g <= 1)){stop("g must be a vector of p-values")}
  
  if(is.null(S)){S <- seq(length(g))}
  if(!is.vector(S) || !is.numeric(S) || !all(is.finite(S))){stop("S must be a vector of finite integers")}
  if(!all(floor(S)==S)){stop("S must be a vector of finite integers")}
  if(!all(S >= 0) || !all(S <= length(g))){stop("S must contain indices between 1 and the total number of variables")}
  S <- unique(S)
  
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  
  type = match.arg(tolower(type), c("fisher", 
                                    "fisher.dep",
                                    "pearson", 
                                    "liptak",
                                    "cauchy", 
                                    "vovk.wang",
                                    "harmonic.dep", 
                                    "harmonic.ind"))
  
  if(type=="fisher" || type=="fisher.dep"){g <- -2*log(g)}
  else if(type=="pearson"){g <- 2*log(1-g)}
  else if(type=="liptak"){g <- qnorm(1-g)}
  else if(type=="cauchy") {g <- tan((0.5-g)*pi)}
  else if (type=="harmonic.dep" || type=="harmonic.ind"){g <- 1/g}
  else if(type=="vovk.wang"){
    if(!is.numeric(r)){stop("r must be a real number")}
    if(r==0){
      g <- -log(g)
    }else if(is.finite(r)){
      g <- -sign(r) * g^r
    }else if (r==Inf){
      g <- max(g)
    }else if (r==-Inf){
      g <- min(g)
    }
  }
  
  if(!all(is.finite(g))){stop("Transformation produced infinite values")}
  
  cvs <- generateCV(length(g), type, alpha, r)
  
  out <- psumTest(g, S, alpha, cvs)
  return(out)
}