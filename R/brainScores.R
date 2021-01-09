#' @title Matrix of Scores for Brain Imaging
#' @description This function computes t-scores for different permutations of brain imaging data.
#' @usage brainScores(copes, mask = NULL, alternative = "two.sided", alpha = 0.05, B = 1000, seed = NULL,
#' truncFrom = 3.2, truncTo = 0, squares = FALSE)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param B number of permutations, including the identity.
#' @param seed seed.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param squares logical, \code{TRUE} to compute squared t-scores.
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{brainScores} returns a list containing the matrix of t-scores \code{G},
#' and the transformed truncation parameters \code{truncFrom} and \code{truncTo}.
#' In \code{G}, columns correspond to variables, and rows to permutations.
#' The first permutation is the identity.
#' @author Anna Vesely.
#' @export


brainScores <- function(copes, mask= NULL, alternative="two.sided", alpha=0.05, B=1000, seed=NULL,
                        truncFrom=3.2, truncTo=0, squares=FALSE){
  
  out <- brainFlip(copes, mask, alternative, alpha, B, seed, truncFrom, truncTo, pvalues=FALSE,
                               type="vovk.wang", r=0, squares, rand=FALSE)
  return(out)
}