#' @title Permutation p-Values for Brain Imaging
#' @description This function computes p-value combinations for different permutations of brain imaging data.
#' A voxel's p-value is calculated by performing the one-sample t test
#' for the null hypothesis that its mean contrast over the different subjects is zero.
#' @usage brainPvals(copes, mask = NULL, alternative = "two.sided", alpha = 0.05, B = 200, 
#'            seed = NULL, truncFrom = alpha, truncTo = max(alpha, 0.5),
#'            type = "vovk.wang", r = 0, rand = FALSE)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain, or character for a Nifti file name.
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
#' \item Cauchy: \code{tan(0.5 - p)/p}
#' \item Vovk and Wang: \code{- sign(r)p^r}
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{brainPvals} returns an object of class \code{sumBrain}, containing
#' \itemize{
#' \item \code{statistics}: numeric matrix of p-values, where columns correspond to voxels inside the brain, and rows to permutations.
#' The first permutation is the identity
#' \item \code{mask}: 3D logical array, where \code{TRUE} values correspond to voxels inside the brain
#' \item \code{alpha}: significance level
#' \item \code{truncFrom}: transformed first truncation parameter
#' \item \code{truncTo}: transformed second truncation parameter
#' }
#' @author Anna Vesely.
#' @examples
#' # use data from the package fMRIdata
#' if(!requireNamespace("fMRIdata", quietly = TRUE)){devtools::install_github("angeella/fMRIdata")}
#' \donttest{
#' library(fMRIdata)
#' data("Auditory_copes")
#' data("Auditory_mask")
#' data("Auditory_clusterTH3_2")
#' 
#' # create object of class sumBrain (combination: Cauchy)
#' res <- brainPvals(copes = Auditory_copes, mask = Auditory_mask, type = "cauchy", seed = 42)
#' 
#' res
#' summary(res)
#' 
#' # confidence bound for the number of true discoveries and the true discovery proportion within clusters
#' # (may require some minutes)
#' out <- clusterAnalysis(res, clusters = Auditory_clusterTH3_2)
#' 
#' # write the TDP map as Nifti file: download mask.nii.gz in the working directory
#' # from https://github.com/angeella/fMRIdata/blob/master/data-raw/AuditoryData
#' library(RNifti)
#' RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz", template = "mask.nii.gz")
#' }
#' @export


brainPvals <- function(copes, mask= NULL, alternative="two.sided", alpha=0.05, B=200, seed=NULL,
                        truncFrom=NULL, truncTo=0.5, type="vovk.wang", r=0, rand=FALSE){
  
  if(!is.null(truncTo) && is.null(truncFrom)){
    truncFrom <- alpha
    truncTo <- max(alpha, truncTo)
  }
  
  out <- brainFlip(copes, mask, alternative, alpha, B, seed, truncFrom, truncTo, pvalues=TRUE,
                   type, r, squares=FALSE, rand)
  return(out)
}



