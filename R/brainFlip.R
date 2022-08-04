#' @title Permutation Statistics for Brain Imaging
#' @description Internal function.
#' It computes test statistics for different permutations of brain imaging data.
#' A voxel's statistic is calculated by performing the one-sample t test
#' for the null hypothesis that its mean contrast over the different subjects is zero.
#' @usage brainFlip(copes, mask, alternative, alpha, B, seed, truncFrom, truncTo, pvalues, type, r, squares, rand)
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
#' @param pvalues logical, \code{TRUE} to use p-values, \code{FALSE} to use t-scores.
#' @param type p-value combination among \code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{vovk.wang} (see details).
#' @param r parameter for Vovk and Wang's p-value transformation.
#' @param squares logical, \code{TRUE} to use squared t-scores.
#' @param rand logical, \code{TRUE} to compute p-values by permutation distribution.
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
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
#' \code{truncTo} should be strictly smaller than 1.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{brainFlip} returns an object of class \code{sumBrain}, containing
#' \itemize{
#' \item \code{statistics}: numeric matrix of statistics, where columns correspond to voxels inside the brain, and rows to permutations.
#' The first permutation is the identity
#' \item \code{mask}: 3D logical array, where \code{TRUE} values correspond to voxels inside the brain
#' \item \code{alpha}: significance level
#' \item \code{truncFrom}: transformed first truncation parameter
#' \item \code{truncTo}: transformed second truncation parameter
#' }
#' @author Anna Vesely.
#' @noRd
#' @importFrom RNifti readNifti
#' @importFrom pARI signTest


brainFlip <- function(copes, mask, alternative, alpha, B, seed, truncFrom, truncTo, pvalues,
                              type, r, squares, rand){
  
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
  
  alternative <- match.arg(tolower(alternative), c("greater", "lower", "two.sided"))
  type <- match.arg(tolower(type), c("fisher", "pearson", "liptak", "edgington", "cauchy", "harmonic", "vovk.wang"))
  
  # check alpha and B
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  if(!is.numeric(B) || !is.finite(B) || B <= 0){stop("B must be a positive integer")}
  B <- ceiling(B)
  if(B < (1/alpha)){stop("1/alpha cannot exceed the number of transformations")}
  
  if(!is.null(seed)){if(!is.numeric(seed) || !is.finite(seed)){stop("seed must be a finite integer")}}
  else{seed <- sample(seq(10^10), 1)}
  set.seed(round(seed))
  
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
  
  scores <- scores[which(mask != 0),]
  if(!is.numeric(scores) || !all(is.finite(scores))){stop("copes should contain numeric values for voxels inside the brain")}
  
  st <- pARI::signTest(X=scores, B=B, alternative=alternative, seed=seed, rand=rand) # sign flipping
  rm(scores)
  
  if(!pvalues){
    G <- rbind(st$Test, t(st$Test_H0))
    option <- ifelse(squares, "squares", alternative)
  }else{
    G <- rbind(st$pv, t(st$pv_H0))
    option <- type
  }
  rm(st)
  
  res <- transf(G, truncFrom, truncTo, option, r)
  out <- sumBrain(res$G, mask, alpha, res$truncFrom, res$truncTo)
  return(out)
}