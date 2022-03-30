#' @title True Discovery Guarantee for Cluster Analysis of Brain Imaging Data
#' @description This function determines a true discovery guarantee for fMRI cluster analysis.
#' It computes confidence bounds for the number of true discoveries and the true discovery proportion
#' within each cluster. The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage brainAnalysis(sumBrain, clusters = NULL, nMax = 50, silent = FALSE)
#' @param sumBrain an object of class sumBrain, as returned by the functions \code{\link{brainScores}} and \code{\link{brainPvals}}.
#' @param clusters 3D numeric array of cluster indices, or character for a Nifti file name.
#' If NULL, the whole brain is considered.
#' @param nMax maximum number of iterations per cluster.
#' @param silent logical, \code{FALSE} to print a summary of active clusters.
#' @return \code{brainAnalysis} returns a list containing \code{summary} (data frame) and
#' \code{TDPmap} (3D numeric array of the true discovery proportions).
#' The data frame \code{summary} contains, for each cluster,
#' \itemize{
#' \item \code{size}: size
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the true discovery proportion
#' \item \code{maxTD}: maximum value of \code{TDP} that could be found under convergence of the algorithm
#' \item \code{dim1}, \code{dim2}, \code{dim3}: coordinates of the center of mass.
#' }
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
#' @references
#' Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.
#' 
#' Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2021). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' Permutation statistics for brain imaging: \code{\link{brainScores}}, \code{\link{brainPvals}}
#' 
#' Suprathreshold clusters: \code{\link{brainClusters}}
#' @importFrom RNifti readNifti
#' @export


brainAnalysis <- function(sumBrain, clusters=NULL, nMax=50, silent=FALSE){
  
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
    clusters <- mask
  }
  
  clusterId <- sort(unique(as.vector(clusters[mask != 0])), decreasing=TRUE) # define number of clusters
  
  M <- matrix(NA, nrow=length(clusterId), ncol=8)
  colnames(M) <- c("size", "TD", "maxTD", "TDP", "maxTDP", "dim1", "dim2", "dim3")
  rownames(M) <- paste("cl", clusterId, sep="")
  
  for(i in seq(length(clusterId))){
    sel <- (clusters == clusterId[i])
    sel[mask==0] <- FALSE
    S <- which(sel[mask != 0])
    out <- sumTest(G, S, alpha, truncFrom, truncTo, nMax)
    
    # cluster summary
    cl <- which(sel, arr.ind=TRUE)
    meanCoord <- colMeans(cl)
    centerInd <- which.min(colSums((t(cl) - meanCoord)^2))
    centerCoord <- cl[centerInd,]
    M[i,] <- c(out$size, out$TD, out$maxTD, round(out$TD/out$size, 3), round(out$maxTD/out$size, 3), centerCoord)
  }
  
  # write TDPmap
  TDPmap <- clusters
  vals <- as.numeric(M[,4])
  vals[length(clusterId)] <- 0
  
  for(i in seq(length(clusterId))){
    TDPmap[clusters==clusterId[i]] <- vals[i]
  }
  
  M <- as.data.frame(M)
  
  if(!silent){
    sel <- M$TD >0
    if(sum(sel) == 0){
      print("No active clusters.")
    }else{
      W <- M[sel,]
      W <- W[order(W$TDP, W$size, decreasing=TRUE),]
      W$converge <- (W$TD == W$maxTD) + 0
      W <- W[,-c(3,5)]
      cat("Active clusters: ", as.character(nrow(W)), ".\n", sep="")
      cat("Clusters with highest TDP:\n")
      cat("\n")
      print(head(W,20))
    }
  }
  
  return(list("summary" = M, "TDPmap" = TDPmap))
}