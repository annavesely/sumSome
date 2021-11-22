#' @title sumBrain Class
#' @description Internal function.
#' It defines a \code{sumBrain} object, storing information on the test statistics
#' for different permutations of brain imaging data.
#' @usage sumBrain(G, mask, alpha, truncFrom, truncTo)
#' @param G numeric matrix of statistics, where columns correspond to voxels inside the brain, and rows to permutations.
#' The first permutation is the identity.
#' @param mask 3D logical array, where \code{TRUE} values correspond to voxels inside the brain.
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @return \code{sumBrain} returns an object of class \code{sumBrain}, containing
#' \code{statistics} (matrix \code{G}), \code{mask}, \code{alpha}, \code{truncFrom} and \code{truncTo}.
#' @author Anna Vesely.
#' @keywords Internal

sumBrain <- function(G, mask, alpha, truncFrom, truncTo){
  val <- list("statistics"=G, "mask"=mask, "alpha"=alpha, "truncFrom"=truncFrom, "truncTo"=truncTo)
  attr(val, "class") <- "sumBrain"
  return(val)
}





#' @export
print.sumBrain <- function(x, ...){
  col <- ncol(x$statistics)
  f <- ifelse(col==1, "1 voxel", paste(as.character(col), " voxels", sep=""))
  cat("A sumBrain object containing brain imaging data with ", f, ".\n", sep="")
  cat("Use clusterAnalysis() for cluster analysis using this object.\n", sep="")
}






#' @export
summary.sumBrain <- function(object, ...){
  alpha <- as.character(object$alpha)
  B <- as.character(nrow(object$statistics))
  
  print(object)
  cat("\n", sep="")
  cat("Significance level: ", alpha, ".\n", sep="")
  cat("Number of permutations: ", B, ".\n", sep="")
}


