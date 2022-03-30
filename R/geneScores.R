#' @title Permutation t-Scores for Gene Expression
#' @description This function computes t-scores for different permutations of gene expression data.
#' A gene's score is calculated by performing the two-sample t test
#' for the null hypothesis that the mean expression value is the same between two populations.
#' @usage geneScores(expr, labels, alternative = "two.sided", alpha = 0.05, B = 200, seed = NULL,
#'            truncFrom = 3.2, truncTo = 0, squares = FALSE)
#' @param expr matrix where rows correspond to genes, and columns to samples.
#' @param labels numeric/character vector with two levels, denoting the population of each sample.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param B number of permutations, including the identity.
#' @param seed seed.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param squares logical, \code{TRUE} to use squared t-scores.
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{geneScores} returns an object of class \code{sumGene}, containing
#' \itemize{
#' \item \code{statistics}: numeric matrix of scores, where columns correspond to genes, and rows to permutations.
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
#' res <- geneScores(expr = expr, labels = labels, alpha = 0.2, seed = 42)
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
#' Vesely, A., Finos, L., and Goeman, J. J. (2021). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' Permutation statistics for gene expression using p-values: \code{\link{genePvals}}
#' 
#' True discovery guarantee for cluster analysis: \code{\link{geneAnalysis}}
#' @export


geneScores <- function(expr, labels, alternative="two.sided", alpha=0.05, B=200, seed=NULL,
                        truncFrom=3.2, truncTo=0, squares=FALSE){
  
  out <- geneFlip(expr, labels, alternative, alpha, B, seed, truncFrom, truncTo, pvalues=FALSE,
                   type="vovk.wang", r=0, squares, rand=FALSE)
  return(out)
}