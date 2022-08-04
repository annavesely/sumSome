#' @title True Discovery Guarantee for p-Value Combinations
#' @description This function determines confidence bounds for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest, when using p-values as test statistics.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage sumPvals(G, S = NULL, alpha = 0.05, truncFrom = NULL, truncTo = 0.5,
#'          type = "vovk.wang", r = 0, nMax = 50)
#' @param G numeric matrix of p-values, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' The first transformation is the identity.
#' @param S vector of indices for the variables of interest (if not specified, all variables).
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, it is set to \code{alpha}.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, p-values are not truncated.
#' @param type p-value combination among \code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{harmonic}, \code{vovk.wang} (see details).
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param nMax maximum number of iterations.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Edgington: \code{p} (Edgington, 1972)
#' \item Fisher: \code{-2log(p)} (Fisher, 1925)
#' \item Pearson: \code{2log(1-p)} (Pearson, 1933)
#' \item Liptak: \code{qnorm(1-p)} (Liptak, 1958; Stouffer et al., 1949)
#' \item Cauchy: \code{tan[(0.5-p)pi]} with \code{pi=3.142} (Liu and Xie, 2020)
#' \item Harmonic mean: \code{1/p} (Wilson, 2019)
#' \item Vovk and Wang: \code{p^r} (\code{log(p)} for \code{r}=0) (Vovk and Wang, 2020)
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1), where
#' \code{B} is the number of data transformations (rows in \code{G}).
#' @return \code{sumPvals} returns an object of class \code{sumObj}, containing
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
#' @references
#' Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.
#' 
#' Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2021). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' True discovery guarantee using generic statistics: \code{\link{sumStats}}
#' 
#' Access a \code{sumObj} object: \code{\link{discoveries}}, \code{\link{tdp}}, \code{\link{fdp}}
#' @export

  
sumPvals <- function(G, S=NULL, alpha=0.05, truncFrom=NULL, truncTo=0.5, type="vovk.wang", r=0, nMax=50){
  
  if(is.null(S)){S <- seq(ncol(G))}
  
  if(is.null(truncFrom)){truncFrom <- alpha}
  
  type = match.arg(tolower(type), c("fisher", "pearson", "liptak", "edgington", "cauchy", "harmonic", "vovk.wang"))
  res <- transf(G, truncFrom, truncTo, type, r)
  rm(G)
  
  out <- sumTest(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)
  return(out)
}