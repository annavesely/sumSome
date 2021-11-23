#' @title Permutation t-Scores for Brain Imaging
#' @description This function computes t-scores for different permutations of brain imaging data.
#' A voxel's score is calculated by performing the one-sample t test
#' for the null hypothesis that its mean contrast over the different subjects is zero.
#' @usage brainScores(copes, mask = NULL, alternative = "two.sided", alpha = 0.05, B = 200,
#'             seed = NULL, truncFrom = 3.2, truncTo = 0, squares = FALSE)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain, or character for a Nifti file name.
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
#' @return \code{brainScores} returns an object of class \code{sumBrain}, containing
#' \itemize{
#' \item \code{statistics}: numeric matrix of t-scores, where columns correspond to voxels inside the brain, and rows to permutations.
#' The first permutation is the identity
#' \item \code{mask}: 3D logical array, where \code{TRUE} values correspond to voxels inside the brain
#' \item \code{alpha}: significance level
#' \item \code{truncFrom}: transformed first truncation parameter
#' \item \code{truncTo}: transformed second truncation parameter
#' }
#' @author Anna Vesely.
#' @examples
#' # # use data from the package fMRIdata
#' # if(!requireNamespace("fMRIdata", quietly = TRUE)){
#' #   devtools::install_github("angeella/fMRIdata")
#' # }
#' # library(fMRIdata)
#' # data("Auditory_copes")
#' # data("Auditory_mask")
#' # data("Auditory_clusterTH3_2")
#' 
#' # # create object of class sumBrain
#' # res <- brainScores(copes = Auditory_copes, mask = Auditory_mask, seed = 42)
#' # res
#' # summary(res)
#' 
#' # # confidence bound for the number of true discoveries and the TDP within clusters
#' # (may require some minutes)
#' # out <- clusterAnalysis(res, clusters = Auditory_clusterTH3_2)
#' 
#' # # write the TDP map as Nifti file: download mask.nii.gz in the working directory
#' # # from https://github.com/angeella/fMRIdata/blob/master/data-raw/AuditoryData
#' # RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz", template = "mask.nii.gz")
#' @references
#' Goeman, J. J., and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science 26 (4) 584-597.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' Permutation statistics for brain imaging using p-values: \code{\link{brainPvals}}
#' 
#' True discovery guarantee for cluster analysis: \code{\link{clusterAnalysis}}
#' 
#' Suprathreshold clusters: \code{\link{findClusters}}
#' @export


brainScores <- function(copes, mask= NULL, alternative="two.sided", alpha=0.05, B=200, seed=NULL,
                        truncFrom=3.2, truncTo=0, squares=FALSE){
  
  out <- brainFlip(copes, mask, alternative, alpha, B, seed, truncFrom, truncTo, pvalues=FALSE,
                               type="vovk.wang", r=0, squares, rand=FALSE)
  return(out)
}