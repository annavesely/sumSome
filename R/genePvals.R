#' @title Permutation p-Values for Gene Expression
#' @description This function computes p-value combinations for different permutations of gene expression data.
#' A gene's p-value is calculated by performing the two-sample t test
#' for the null hypothesis that the mean expression value is the same between two populations.
#' @usage genePvals(expr, labels, alternative = "two.sided", alpha = 0.05, B = 200, seed = NULL,
#'           truncFrom = NULL, truncTo = 0.5, type = "vovk.wang", r = 0, rand = FALSE)
#' @param expr matrix where rows correspond to genes, and columns to samples.
#' @param labels numeric/character vector with two levels, denoting the population of each sample.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param B number of permutations, including the identity.
#' @param seed seed.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, it is set to \code{alpha}.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, p-values are not truncated.
#' @param type p-value combination among \code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{vovk.wang} (see details).
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param rand logical, \code{TRUE} to compute p-values by permutation distribution.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Edgington: \code{-p}
#' \item Fisher: \code{-log(p)}
#' \item Pearson: \code{log(1-p)}
#' \item Liptak: \code{-qnorm(p)}
#' \item Cauchy: \code{tan[(0.5 - p)pi]} with \code{pi}=3.142
#' \item Vovk and Wang: \code{- sign(r)p^r}
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{genePvals} returns an object of class \code{sumGene}, containing
#' \itemize{
#' \item \code{statistics}: numeric matrix of p-values, where columns correspond to genes, and rows to permutations.
#' The first permutation is the identity
#' \item \code{alpha}: significance level
#' \item \code{truncFrom}: transformed first truncation parameter
#' \item \code{truncTo}: transformed second truncation parameter
#' }
#' @author Anna Vesely.
#' @examples
#' # simulate 20 samples of 100 genes
#' set.seed(42)
#' expr <- matrix(c(rnorm(1000, mean = 0, sd = 10), rnorm(1000, mean = 13, sd = 10)), ncol = 20)
#' rownames(expr) <- seq(100)
#' labels <- rep(c(1,2), each = 10)
#' 
#' # simulate pathways
#' pathways <- lapply(seq(3), FUN = function(x) sample(rownames(expr), 3*x))
#' 
#' # create object of class sumGene
#' res <- genePvals(expr = expr, labels = labels, alpha = 0.2, seed = 42, type = "liptak")
#' res
#' summary(res)
#' 
#' # confidence bound for the number of true discoveries and the TDP within pathways
#' out <- geneAnalysis(res, pathways = pathways)
#' out
#' @references
#' Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.
#' 
#' Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' Permutation statistics for gene expression using t scores: \code{\link{geneScores}}
#' 
#' True discovery guarantee for cluster analysis: \code{\link{geneAnalysis}}
#' @export


genePvals <- function(expr, labels, alternative="two.sided", alpha=0.05, B=200, seed=NULL,
                       truncFrom=NULL, truncTo=0.5, type="vovk.wang", r=0, rand=FALSE){
  
  if(is.null(truncFrom)){truncFrom <- alpha}
  
  out <- geneFlip(expr, labels, alternative, alpha, B, seed, truncFrom, truncTo, pvalues=TRUE,
                   type, r, squares=FALSE, rand)
  return(out)
}



