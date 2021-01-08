#' @title True Discovery Guarantee for Brain Imaging
#' @description Internal function, called in \code{sumBrain} and \code{sumBrain.pvals}.
#' It determines a true discovery guarantee for fMRI cluster analysis.
#' @usage sumBrain.internal(copes, mask, clusters, thr, alternative, alpha, B, seed, truncFrom, truncTo,
#'                   pvalues, type, r, squares, rand, nMax, silent)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain, or character for a Nifti file name.
#' @param clusters 3D numeric array of cluster indices, or character for a Nifti file name.
#' @param thr threshold used to compute the clusters: adjacent t-scores more extreme than \code{thr}
#' are combined into a single cluster.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param alpha significance level.
#' @param B number of permutations, including the identity.
#' @param seed seed.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param pvalues logical, \code{TRUE} to use p-values, \code{FALSE} to use t-scores.
#' @param type transformation of p-values (\code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{vovk.wang})
#' @param r parameter for Vovk and Wang's p-value transformation.
#' @param squares logical, \code{TRUE} to use squared t-scores.
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
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{sumBrain.internal} returns a list containing \code{summary} (matrix),
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
#' @importFrom RNifti readNifti
#' @importFrom pARI signTest
#' @importFrom ARIbrain cluster_threshold


brainClusters <- function(copes, mask, thr, alternative){
  
  # check copes
  if(!is.list(copes)){stop("copes should be a list of arrays")}
  n <- length(copes)
  if(n==0){stop("copes should be a list of arrays")}
  imgDim <- dim(copes[[1]]) #(91,109,91)
  
  # check mask
  if(!is.null(mask)){
    if(!is.character(mask) && !is.array(mask)){stop("mask must be an array or a path")}
    if(is.character(mask)){mask = RNifti::readNifti(mask)}
    if(!all(dim(mask) == imgDim)){stop("incompatible dimensions of mask and copes")}
  }else{
    mask <- array(1, imgDim)
  }
  
  if(is.null(thr)){
    if(!is.null(mask)){clusters <- array(mask, dim(mask))}
    else{stop("Please insert mask, threshold value or cluster map")}
  }
  
  # check threshold
  if(!is.null(thr) && !(is.numeric(thr) && is.finite(thr))){stop("thr must be a finite number")}
  
  alternative <- match.arg(tolower(alternative), c("greater", "lower", "two.sided"))
  
  # create image
  img <- array(NA, c(imgDim, n))
  for (i in seq(n)) {
    if(!(all(dim(copes[[i]]) == imgDim))){stop("incompatible copes dimensions")}
    img[,,,i] <- copes[[i]]
  }
  rm(copes)
  
  # matrix of data (rows = variables, columns = observations)
  scores <- matrix(img, nrow=(imgDim[1] * imgDim[2] * imgDim[3]), ncol=n)
  scores[mask==0,] <- NA
  rm(img)
  
  # if needed, create clusters
  if(is.null(clusters)){
    # maps of t-statistics
    tMap <- pARI::signTest(scores, 1, alternative, rand = FALSE)$Test
    tMap <- array(tMap, imgDim)
    tMap[mask==0] <- 0
    if(alternative == "greater"){clusters <- ARIbrain::cluster_threshold(tMap > thr)}
    else if(alternative == "two.sided"){clusters <- ARIbrain::cluster_threshold(abs(tMap) > abs(thr))}
    else{clusters <- ARIbrain::cluster_threshold(-tMap < thr)}
  }
  
  return(clusters)
}