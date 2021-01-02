#' @title True Discovery Guarantee (General Case)
#' @description Internal function, called in \code{sumSome}, \code{sumSome.pvals} and \code{sumBrain.internal}.
#' It determines a lower confidence bound for the number of true discoveries within a set of interest.
#' The bound remains valid under post-hoc selection.
#' @usage sum.internal(G, S, alpha, truncFrom, truncTo, nMax)
#' @param G numeric matrix of statistics, where columns correspond to variables, and rows to data transformations (e.g. permutations).
#' Extreme values are the greatest.
#' @param S vector of indices for the variables of interest.
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @param nMax maximum number of iterations.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @return \code{sum.internal} returns a list containing \code{summary} (vector) and \code{iterations} (number of iterations).
#' \code{summary} contains:
#' \itemize{
#' \item \code{size}: size of \code{S}
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the true discovery proportion in \code{S}
#' \item \code{maxTD}: maximum value of \code{TDP} that could be found under convergence of the algorithm
#' }
#' @author Anna Vesely.
#' @keywords internal


sum.internal <- function(G, S, alpha, truncFrom, truncTo, nMax){
  
  B <- nrow(G)
  f <- ncol(G)
  
  if(!is.vector(S) || !is.numeric(S)){stop("S must be a vector of finite integers")}
  if(!all(floor(S)==S)){stop("S must be a vector of finite integers")}
  if(!all(S >= 0) || !all(S <= f)){stop("S must contain indices between 1 and the total number of variables")}
  S <- unique(S)
  s <- length(S)
  
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  if(B < (1/alpha)){stop("1/alpha cannot exceed the number of transformations")}
  k <- ceiling((1-alpha)*B)
  
  if(!is.numeric(nMax) || !is.finite(nMax)){stop("nMax must be a finite number")}
  if(nMax < 0){stop("nMax must be a finite number")}
  nMax <- max(round(nMax), 0)
  
  G <- cbind(G[,S], G[,-S])
  S <- seq(s)
  truncation <- (!is.null(truncFrom) && !is.null(truncTo))
  
  if(truncation){
    # collapse all columns having obs=truncTo (which do not improve the bounds)
    collapse <- which(G[1,] == truncTo)
    collapse <- collapse[collapse > s]
    if(length(collapse) > 0){G <- as.matrix(cbind(G[,-collapse], rowSums(as.matrix(G[,collapse]))))}
    
    # remove columns where all values (except eventually the first one) are = truncTo (that do not worsen the bounds)
    rem <- which(apply(G, 2, function(x) permMin(x, B, truncTo)))
    rem <- rem[rem > s]
    if(length(rem) > 0){G <- as.matrix(G[,-rem])}
    f <- ncol(G)
    
    # reorder indices in increasing order (obs=truncTo are added at the beginning)
    minObs <- which(G[1,] == truncTo)
    otherObs <- setdiff(seq(f), minObs)
    otherObs <- otherObs[order(G[1, otherObs], decreasing=FALSE)]
    indices <- c(minObs, otherObs)
  }else{
    indices <- order(G[1,], decreasing=FALSE)
  }
  
  # centered test statistics with observations in increasing order
  D <- G[,indices]
  D <- sweep(D, 2, D[1,])
  
  # centered test statistics where each row is in ascending order
  o <- t(apply(D, 1, order, decreasing=TRUE))
  R <- t(sapply(seq(B), function(x) D[x, o[x,]]))
  
  # matrix of indices in R
  I <- matrix(rep(indices, B), ncol=f, byrow=TRUE)
  I <- t(sapply(seq(B), function(x) I[x, o[x,]]))
  f <- ncol(I)
  
  b <- bisectionTD(D, I, R, s, f, k, B, nMax)
  summary <- c(s, b$TDmin, b$TDmax, round(b$TDmin * 100/s, 1), round(b$TDmax * 100/s, 1))
  names(summary) <- c("size", "TD", "maxTD", "TDP", "maxTDP")
  out <- list("summary" = summary, "iterations" = b$BAB)
  return(out)
}