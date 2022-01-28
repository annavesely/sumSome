#' @title True Discovery Guarantee for Cluster Analysis of Gene Expression Data
#' @description This function determines a true discovery guarantee for gene cluster analysis.
#' It computes confidence bounds for the number of true discoveries and the true discovery proportion
#' within each cluster. The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage geneAnalysis(sumGene, clusters = NULL, nMax = 50, silent = FALSE)
#' @param sumGene an object of class sumGene, as returned by the functions \code{\link{geneScores}} and \code{\link{genePvals}}.
#' @param pathways list of character vectors containing gene names (one vector per pathway). If NULL, the whole gene set is considered.
#' @param nMax maximum number of iterations per cluster.
#' @param silent logical, \code{FALSE} to print a summary of active pathways.
#' @return \code{geneAnalysis} returns a list matrix containing, for each pathway,
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
#' expr <- matrix(c(rnorm(1000, mean = 0, sd = 10), rnorm(1000, mean = 13, sd = 10)), ncol = 20)
#' rownames(expr) <- seq(100)
#' labels <- rep(c(1,2), each = 10)
#' 
#' # simulate pathways
#' pathways <- lapply(seq(3), FUN = function(x) sample(rownames(expr), 3*x))
#' 
#' # create object of class sumGene
#' res <- geneScores(expr = expr, labels = labels, alpha = 0.2, seed = 42)
#' res
#' summary(res)
#' 
#' # confidence bound for the number of true discoveries and the TDP within pathways
#' out <- geneAnalysis(res, pathways = pathways)
#' out
#' @references
#' Goeman, J. J. and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science, 26(4):584-597.
#' 
#' Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
#' 
#' Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
#' @seealso
#' Permutation statistics for gene expression: \code{\link{geneScores}}, \code{\link{genePvals}}
#' @export


geneAnalysis <- function(sumGene, pathways=NULL, nMax=50, silent=FALSE){
  
  if(class(sumGene) != "sumGene"){stop("sumGene should be an object of class sumGene")}
  
  G <- sumGene$statistics
  alpha <- sumGene$alpha
  truncFrom <- sumGene$truncFrom
  truncTo <- sumGene$truncTo
  rm(sumGene)
  genes <- colnames(G)
  
  # check clusters (if clusters = NULL and thershold != NULL, clusters are computed later)
  if(!is.null(pathways)){
    if(!is.list(pathways) || length(pathways)==0){stop("pathways must be a list")}
  }else{
    pathways <- list("all"=genes)
  }
  
  if(length(names(pathways))==0){names(pathways) <- paste0("pathway", seq(length(pathways)))}
  pathId <- names(pathways)
  
  M <- matrix(NA, nrow=length(pathId), ncol=5)
  colnames(M) <- c("size", "TD", "maxTD", "TDP", "maxTDP")
  rownames(M) <- pathId
  
  for(i in seq(length(pathId))){
    sel <- sapply(genes, FUN=function(x) x %in% pathways[[i]])
    S <- which(sel)
    out <- sumTest(G, S, alpha, truncFrom, truncTo, nMax)
    
    # pathway summary
    M[i,] <- c(out$size, out$TD, out$maxTD, round(out$TD/out$size, 3), round(out$maxTD/out$size, 3))
  }
  
  if(!silent){
    sel <- M[,2]>0
    if(sum(sel) == 0){
      print("No active pathways")
    }else{
      W <- M[sel,,drop=FALSE]
      W <- W[order(W[,1], W[,4], decreasing=TRUE),,drop=FALSE]
      converge <- (W[,2] == W[,3]) + 0
      W <- cbind(W[,-c(3,5),drop=FALSE], converge)
      print(W)
    }
  }
  
  return(M)
}