#' @title True Discovery Guarantee for Brain Imaging
#' @description This function determines a true discovery guarantee for fMRI cluster analysis, using t-scores.
#' @usage sumSomeBrain(copes, mask = NULL, clusters = NULL, thr = 3.2, alternative = "two.sided", alpha = 0.05,
#'              B = 1000, seed = NULL, truncFrom = thr, truncTo = 0, squares = FALSE, nMax = 10000,
#'              silent = FALSE)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain.
#' @param clusters 3D numeric array of cluster indices, or character for a Nifti file name.
#' @param thr threshold used to compute the clusters: adjacent t-scores more extreme than \code{thr}
#' are combined into a single cluster.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param B number of permutations.
#' @param seed seed.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param squares logical, \code{TRUE} to use squared t-scores.
#' @param nMax maximum number of iterations.
#' @param silent logical, \code{FALSE} to print the summary.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' @return \code{sumSomeBrain} returns a list containing \code{summary} (matrix),
#' \code{clusters} (3D numeric array of cluster indices), and
#' \code{TDPmap} (3D numeric array of the true discovery proportions).
#' The matrix \code{summary} contains, for each cluster:
#' \itemize{
#' \item \code{size}: size of the cluster
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the true discovery proportion
#' \item \code{maxTD}: maximum value of \code{TDP} that could be found under convergence of the algorithm
#' \item \code{dim1}, \code{dim2}, \code{dim3}: coordinates of the center of mass.
#' }
#' @author Anna Vesely.
#' @examples
#' devtools::install_github("angeella/fMRIdata")
#' library(fMRIdata)
#' library(RNifti)
#' 
#' out <- sumSomeBrain(copes = Auditory_copes, mask = Auditory_mask, clusters = Auditory_clusterTH3_2,
#'                     seed = 42, nMax = 30)
#' 
#' RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz")
#' @export


sumSomeBrain <- function(copes, mask=NULL, clusters=NULL, thr=3.2, alternative="two.sided",
                         alpha=0.05, B=1000, seed=NULL, truncFrom=thr, truncTo=0, squares=FALSE,
                         nMax=10000, silent=FALSE){
  
  out <- sumSomeBrain.internal(copes, mask, clusters, thr, alternative,
                               alpha, B, seed, truncFrom, truncTo, pvalues=FALSE,
                               type="vovk.wang", r=1, squares, rand=FALSE, nMax, silent)
  
  return(out)
}




