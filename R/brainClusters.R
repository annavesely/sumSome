#' @title Suprathreshold Clusters for Brain Imaging
#' @description This function determines spatially connected clusters, where t-scores are more extreme than a given threshold.
#' @usage brainClusters(copes, mask = NULL, thr = 3.2, alternative = "two.sided", silent = FALSE)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain, or character for a Nifti file name.
#' @param thr threshold.
#' @param alternative direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}).
#' @param silent logical, \code{FALSE} to print the number of clusters.
#' @return \code{brainClusters} returns a 3D numeric array, with integer values corresponding to clusters,
#' and 0 to other voxels.
#' @author Anna Vesely.
#' @examples
#' # simulate 20 copes with dimensions 10x10x10
#' set.seed(42)
#' copes <- list()
#' for(i in seq(20)){copes[[i]] <- array(rnorm(10^3, mean = -10, sd = 30), dim=c(10,10,10))}
#' 
#' # cluster map where t scores are grater than 2.8, in absolute value
#' thr <- 2.8
#' cl <- brainClusters(copes = copes, thr = thr)
#' 
#' # create object of class sumBrain
#' res <- brainScores(copes = copes, alpha = 0.2, seed = 42, truncFrom = thr)
#' res
#' summary(res)
#' 
#' # confidence bound for the number of true discoveries and the TDP within clusters
#' out <- brainAnalysis(res, clusters = cl$clusters)
#' out$summary
#' @seealso
#' Permutation statistics for brain imaging: \code{\link{brainScores}}, \code{\link{brainPvals}}
#' 
#' True discovery guarantee for cluster analysis: \code{\link{brainAnalysis}}
#' @export
#' @importFrom RNifti readNifti
#' @importFrom pARI signTest
#' @importFrom ARIbrain cluster_threshold


brainClusters <- function(copes, mask=NULL, thr=3.2, alternative="two.sided", silent=FALSE){
  
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
  
  # if threshold is NULL, clusters = mask
  if(is.null(thr)){clusters <- array(mask, dim(mask))}else{clusters <- NULL}
  
  # check threshold
  if(!is.null(thr) && !(is.numeric(thr) && is.finite(thr))){stop("thr must be a finite number")}
  
  alternative <- match.arg(tolower(alternative), c("greater", "lower", "two.sided"))
  
  # matrix of data (rows = variables, columns = observations)
  scores <- matrix(NA, nrow=(imgDim[1] * imgDim[2] * imgDim[3]), ncol=n)
  
  for (i in seq(n)) {
    if(!(all(dim(copes[[i]]) == imgDim))){stop("incompatible copes dimensions")}
    scores[,i] <- as.vector(copes[[i]])
  }
  rm(copes)
  
  # NA for voxels not in the brain
  scores[mask==0,] <- NA
  
  # if needed, create clusters
  if(is.null(clusters)){
    # maps of t-statistics
    tMap <- pARI::signTest(scores, 1, alternative, rand = FALSE)$Test
    tMap <- array(tMap, imgDim)
    tMap[mask==0] <- 0
    
    if(alternative == "two.sided"){
      tMap <- abs(tMap)
      thr <- abs(thr)
    }
    else if(alternative == "lower"){
      tMap <- -tMap
      thr <- -thr
    }
    
    if(max(tMap) <= thr){clusters <- array(0, imgDim)}else{clusters <- ARIbrain::cluster_threshold(tMap > thr)}
  }
  
  ncl <- max(as.vector(clusters[mask != 0]))
  out <- list("nClusters"=ncl, "clusters"=clusters)
  
  if(!silent){
    cat("Number of clusters: ", as.character(ncl), ".\n", sep="")
  }
  
  return(out)
}