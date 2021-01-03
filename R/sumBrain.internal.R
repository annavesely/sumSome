#' @title True Discovery Guarantee for Brain Imaging
#' @description Internal function, called in \code{sumBrain} and \code{sumBrain.pvals}.
#' It determines a true discovery guarantee for fMRI cluster analysis.
#' @usage sumBrain.internal(copes, mask, clusters, thr, alternative, alpha, B, seed, truncFrom, truncTo,
#'                   pvalues, type, r, squares, rand, nMax, silent)
#' @param copes list of 3D numeric arrays (contrasts maps for each subject).
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain.
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
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @details Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Edgington: \code{-p}
#' \item Fisher: \code{-log(p)}
#' \item Pearson: \code{log(1-p)}
#' \item Liptak: \code{-qnorm(p)}
#' \item Cauchy: \code{tan(0.5 - p)/p}
#' \item Vovk and Wang: \code{- sign(r)p^r}
#' }
#' @details Pearson's and Liptak's transformations produce infinite values in \code{1}.
#' For such transformations, \code{truncTo} is coerced to be not greater than \code{1 -  .Machine$double.eps}.
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


sumBrain.internal <- function(copes, mask, clusters, thr, alternative,
                                  alpha, B, seed, truncFrom, truncTo, pvalues,
                                  type, r, squares, rand, nMax, silent){
  
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
  
  # check clusters (if clusters = NULL and thershold != NULL, clusters are computed later)
  if(!is.null(clusters)){
    if(!is.character(clusters) && !is.array(clusters)){stop("clusters must be an array or a path")}
    if(is.character(clusters)){clusters = RNifti::readNifti(clusters)}
  }else if(is.null(thr)){
    if(!is.null(mask)){clusters <- array(mask, dim(mask))}
    else{stop("Please insert mask, threshold value or cluster map")}
  }
  
  # check threshold
  if(!is.null(thr) && !(is.numeric(thr) && is.finite(thr))){stop("thr must be a finite number")}
  
  alternative <- match.arg(tolower(alternative), c("greater", "lower", "two.sided"))
  type <- match.arg(tolower(type), c("fisher", "pearson", "liptak", "edgington", "cauchy", "vovk.wang"))
  
  # check alpha and B
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  if(!is.numeric(B) || !is.finite(B) || B <= 0){stop("B must be a positive integer")}
  B <- round(B)
  if(B < (1/alpha)){stop("1/alpha cannot exceed the number of transformations")}
  
  if(!is.null(seed)){
    if(!is.numeric(seed) || !is.finite(seed)){stop("seed must be a finite integer")}
    set.seed(round(seed))
  }
  
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
    if(alternative == "greater"){clusters <- pARI::cluster_threshold(tMap > thr)}
    else if(alternative == "two.sided"){clusters <- pARI::cluster_threshold(abs(tMap) > abs(thr))}
    else{clusters <- pARI::cluster_threshold(-tMap < thr)}
    
    rm(tMap)
  }
  
  scores <- scores[which(mask != 0),]
  st <- pARI::signTest(scores, B, alternative, rand) # sign flipping
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
  rm(G)
  
  clusterId <- sort(unique(as.vector(clusters[mask != 0])), decreasing=TRUE) # define number of clusters
  M <- matrix(NA, nrow=length(clusterId), ncol=8)
  colnames(M) <- c("size", "TD", "maxTD", "TDP", "maxTDP", "dim1", "dim2", "dim3")
  rownames(M) <- clusterId
  TDPmap <- array(0, imgDim)
  
  for(i in seq(length(clusterId))){
    sel <- (clusters == clusterId[i])
    sel[mask==0] <- FALSE
    S <- which(sel[mask != 0])
    
    out <- sum.internal(res$G, S, alpha, res$truncFrom, res$truncTo, nMax)$summary
    TDPmap <- TDPmap + (sel * round(out["TD"] * 100/ out["size"]))
    
    # cluster summary
    cl <- which(sel, arr.ind=TRUE)
    meanCoord <- colMeans(cl)
    centerInd <- which.min(colSums((t(cl) - meanCoord)^2))
    centerCoord <- cl[centerInd,]
    M[i,] <- c(out, centerCoord)
  }
  
  if(!silent){print(M)}
  return(list("summary" = M, "clusters" = clusters, "TDPmap" = TDPmap))
}