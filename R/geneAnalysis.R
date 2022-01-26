#' @title True Discovery Guarantee for Cluster Analysis of Gene Expression Data
#' @description This function determines a true discovery guarantee for gene cluster analysis.
#' It computes confidence bounds for the number of true discoveries and the true discovery proportion
#' within each cluster. The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage geneAnalysis(sumGene, clusters, nMax = 50, silent = FALSE)
#' @param sumGene an object of class sumGene, as returned by the functions \code{\link{geneScores}} and \code{\link{genePvals}}.
#' @param clusters numeric vector of cluster indices. If NULL, the whole gene set is considered.
#' @param nMax maximum number of iterations per cluster.
#' @param silent logical, \code{FALSE} to print the summary.
#' @return \code{geneAnalysis} returns a list containing \code{summary} (matrix) and
#' \code{TDPmap} (numeric vector of the true discovery proportions).
#' The matrix \code{summary} contains, for each cluster,
#' \itemize{
#' \item \code{size}: size
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the true discovery proportion
#' \item \code{maxTD}: maximum value of \code{TDP} that could be found under convergence of the algorithm.
#' }
#' @author Anna Vesely.
#' @examples
#' # simulate 20 samples of 100 genes
#' set.seed(42)
#' expr <- matrix(c(rnorm(1000, mean = 0, sd = 10), rnorm(1000, mean = 10, sd = 10)), ncol = 20)
#' labels <- rep(c(1,2), each = 10)
#' 
#' # simulate clusters
#' clusters <- sample(x = seq(5), size = 100, replace = TRUE)
#' 
#' # create object of class sumGene
#' res <- geneScores(expr = expr, labels = labels, alpha = 0.2, seed = 42)
#' res
#' summary(res)
#' 
#' # confidence bound for the number of true discoveries and the TDP within clusters
#' out <- geneAnalysis(res, clusters = clusters)
#' @references
#' Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.
#' 
#' Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' Permutation statistics for gene expression: \code{\link{geneScores}}, \code{\link{genePvals}}
#' @export


geneAnalysis <- function(sumGene, clusters, nMax=50, silent=FALSE){
  
  if(class(sumGene) != "sumGene"){stop("sumGene should be an object of class sumGene")}
  
  G <- sumGene$statistics
  alpha <- sumGene$alpha
  truncFrom <- sumGene$truncFrom
  truncTo <- sumGene$truncTo
  rm(sumGene)
  
  # check clusters (if clusters = NULL and thershold != NULL, clusters are computed later)
  if(!is.null(clusters)){
    if(!is.vector(clusters)){stop("clusters must be a vector")}
    if(length(clusters) != ncol(G)){stop("Incorrect clusters dimensions")}
  }else{
    clusters <- rep(1, ncol(G))
  }
  
  clusterId <- sort(unique(clusters), decreasing=TRUE) # define number of clusters
  
  M <- matrix(NA, nrow=length(clusterId), ncol=5)
  colnames(M) <- c("size", "TD", "maxTD", "TDP", "maxTDP")
  rownames(M) <- paste("cl", clusterId, sep="")
  
  for(i in seq(length(clusterId))){
    S <- which(clusters == clusterId[i])
    out <- sumTest(G, S, alpha, truncFrom, truncTo, nMax)
    
    # cluster summary
    M[i,] <- c(out$size, out$TD, out$maxTD, round(out$TD/out$size, 3), round(out$maxTD/out$size, 3))
  }
  
  # write TDPmap
  TDPmap <- clusters
  vals <- as.numeric(M[,4])
  
  for(i in seq(length(clusterId))){
    TDPmap[clusters==clusterId[i]] <- vals[i]
  }
  
  if(!silent){print(M)}
  return(list("summary" = M, "TDPmap" = TDPmap))
}