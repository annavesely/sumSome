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
#' @export


brainDiscoveries <- function(sumBrain, clusters, nMax=1000, silent=FALSE){
  
  if(class(sumBrain) != "sumBrain"){stop("sumBrain should be an object of class sumBrain")}
  
  G <- sumBrain$statistics
  mask <- sumBrain$mask
  alpha <- sumBrain$alpha
  truncFrom <- sumBrain$truncFrom
  truncTo <- sumBrain$truncTo
  rm(sumBrain)
  
  imgDim <- dim(mask)
  
  # check clusters (if clusters = NULL and thershold != NULL, clusters are computed later)
  if(!is.null(clusters)){
    if(!is.character(clusters) && !is.array(clusters)){stop("clusters must be an array or a path")}
    if(is.character(clusters)){clusters = RNifti::readNifti(clusters)}
    if(!all(dim(clusters) == imgDim)){stop("Incorrect clusters dimensions")}
  }else{
    clusters <- array(mask, imgDim)
  }
  
  clusterId <- sort(unique(as.vector(clusters[mask != 0])), decreasing=TRUE) # define number of clusters
  
  M <- matrix(NA, nrow=length(clusterId), ncol=8)
  colnames(M) <- c("size", "TD", "maxTD", "TDP", "maxTDP", "dim1", "dim2", "dim3")
  rownames(M) <- clusterId
  TDPmap <- array(0, imgDim)
  
  for(i in seq(length(clusterId))){
    sel <- (clusters == clusterId[i])
    sel[mask==0] <- FALSE
    S <- which(sel[mask != 0])
    
    out <- sumTest(G, S, alpha, truncFrom, truncTo, nMax)
    TDPmap <- TDPmap + (sel * (out$TDP * 100))
    
    # cluster summary
    cl <- which(sel, arr.ind=TRUE)
    meanCoord <- colMeans(cl)
    centerInd <- which.min(colSums((t(cl) - meanCoord)^2))
    centerCoord <- cl[centerInd,]
    M[i,] <- c(out$size, out$TD, out$maxTD, out$TDP, out$maxTDP, centerCoord)
  }
  
  if(!silent){print(M)}
  return(list("summary" = M, "TDPmap" = TDPmap))
}