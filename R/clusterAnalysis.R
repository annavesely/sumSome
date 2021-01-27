#' @title True Discovery Guarantee for Cluster Analysis
#' @description This function determines a true discovery guarantee for fMRI cluster analysis.
#' It computes confidence bounds for the number of true discoveries and the true discovery proportion
#' within each cluster. The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage clusterAnalysis(sumBrain, clusters, nMax = 1000, silent = FALSE)
#' @param sumBrain an object of class sumBrain, as returned by the functions \code{\link{brainScores}} and \code{\link{brainPvals}}.
#' @param clusters 3D numeric array of cluster indices, or character for a Nifti file name.
#' If NULL, the whole brain is considered.
#' @param nMax maximum number of iterations per cluster.
#' @param silent logical, \code{FALSE} to print the summary.
#' @return \code{clusterAnalysis} returns a list containing \code{summary} (matrix) and
#' \code{TDPmap} (3D numeric array of the true discovery proportions).
#' The matrix \code{summary} contains, for each cluster,
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
#' # if needed, install the package fMRIdata from Github
#' devtools::install_github("angeella/fMRIdata")
#' 
#' library(fMRIdata)
#' data("Auditory_copes")
#' data("Auditory_mask")
#' data("Auditory_clusterTH3_2")
#' 
#' # create object of class sumBrain
#' res <- brainScores(copes = Auditory_copes, mask = Auditory_mask, B = 200, seed = 42)
#' 
#' res
#' summary(res)
#' 
#' # confidence bound for the number of true discoveries and the true discovery proportion within clusters
#' # (may require some minutes)
#' out <- clusterAnalysis(res, clusters = Auditory_clusterTH3_2, nMax=50, silent=FALSE)
#' 
#' # write the TDP map as Nifti file
#' library(RNifti)
#' RNifti::writeNifti(out$TDPmap, file = "TDPmap.nii.gz")
#' @importFrom RNifti readNifti
#' @export


clusterAnalysis <- function(sumBrain, clusters, nMax=1000, silent=FALSE){
  
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
    TDPmap <- TDPmap + (sel * round(out$TD * 100/ out$size))
    
    # cluster summary
    cl <- which(sel, arr.ind=TRUE)
    meanCoord <- colMeans(cl)
    centerInd <- which.min(colSums((t(cl) - meanCoord)^2))
    centerCoord <- cl[centerInd,]
    M[i,] <- c(out$size, out$TD, out$maxTD, round(out$TD/out$size, 3), round(out$maxTD/out$size, 3), centerCoord)
  }
  
  if(!silent){print(M)}
  return(list("summary" = M, "TDPmap" = TDPmap))
}