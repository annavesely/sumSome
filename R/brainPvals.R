#' @title True Discovery Guarantee for Brain Imaging with p-Values
#' @description This function determines a true discovery guarantee for fMRI cluster analysis, using p-values.
#' @usage sumBrain.pvals(copes, mask = NULL, clusters = NULL, thr = 3.2, alternative = "two.sided",
#'                alpha = 0.05, B = 1000, seed = NULL, truncFrom = alpha, truncTo = 1, type = "fisher",
#'                r = 1, rand = FALSE, nMax = 10000, silent = FALSE)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain.
#' @param clusters 3D numeric array of cluster indices, or character for a Nifti file name.
#' @param thr threshold used to compute the clusters: adjacent t-scores more extreme than \code{thr}
#' are combined into a single cluster.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param B number of permutations, including the identity.
#' @param seed seed.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, p-values are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, p-values are not truncated.
#' @param type transformation of p-values (\code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{vovk.wang})
#' @param r parameter for Vovk and Wang's p-value transformation.
#' @param rand logical, \code{TRUE} to compute p-values by permutation distribution.
#' @param nMax maximum number of iterations.
#' @param silent logical, \code{FALSE} to print the summary.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Edgington: \code{-p}
#' \item Fisher: \code{-log(p)}
#' \item Pearson: \code{log(1-p)}
#' \item Liptak: \code{-qnorm(p)}
#' \item Cauchy: \code{tan(0.5 - p)/p}
#' \item Vovk and Wang: \code{- sign(r)p^r}
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{sumBrain.pvals} returns a list containing \code{summary} (matrix),
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
#' # if needed, install the package fMRIdata from Github
#' devtools::install_github("angeella/fMRIdata")
#' 
#' library(fMRIdata)
#' data("Auditory_copes")
#' data("Auditory_mask")
#' data("Auditory_clusterTH3_2")
#' 
#' # the following requires some minutes
#' out <- sumBrain.pvals(copes = Auditory_copes, mask = Auditory_mask, clusters = Auditory_clusterTH3_2,
#'                       B = 100, seed = 42, type = "fisher", nMax = 30)
#' 
#' # write the TDP map as Nifti file
#' library(RNifti)
#' RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz")
#' @export


brainPvals <- function(copes, mask= NULL, alternative="two.sided", alpha=0.05, B=1000, seed=NULL,
                        truncFrom=alpha, truncTo=NULL, type="vovk.wang", r=0, rand=FALSE){
  
  out <- brainFlip(copes, mask, alternative, alpha, B, seed, truncFrom, truncTo, pvalues=TRUE,
                   type, r, squares=FALSE, rand)
  return(out)
}



