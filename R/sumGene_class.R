#' @title sumGene Class
#' @description Internal function.
#' It defines a \code{sumGene} object, storing information on the test statistics
#' for different permutations of gene expression data.
#' @usage sumGene(G, alpha, truncFrom, truncTo)
#' @param G numeric matrix of statistics, where columns correspond to genes, and rows to permutations.
#' The first permutation is the identity.
#' @param alpha significance level.
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated.
#' If \code{NULL}, statistics are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, statistics are not truncated.
#' @return \code{sumGene} returns an object of class \code{sumGene}, containing
#' \code{statistics} (matrix \code{G}), \code{alpha}, \code{truncFrom} and \code{truncTo}.
#' @author Anna Vesely.
#' @noRd

sumGene <- function(G, alpha, truncFrom, truncTo){
  val <- list("statistics"=G, "alpha"=alpha, "truncFrom"=truncFrom, "truncTo"=truncTo)
  attr(val, "class") <- "sumGene"
  return(val)
}





#' @export
print.sumGene <- function(x, ...){
  col <- ncol(x$statistics)
  f <- ifelse(col==1, "1 gene", paste(as.character(col), " genes", sep=""))
  cat("A sumGene object containing gene expression data with ", f, ".\n", sep="")
  cat("Use geneAnalysis() for cluster analysis using this object.\n", sep="")
}






#' @export
summary.sumGene <- function(object, ...){
  alpha <- as.character(object$alpha)
  B <- as.character(nrow(object$statistics))
  
  print(object)
  cat("\n", sep="")
  cat("Significance level: ", alpha, ".\n", sep="")
  cat("Number of permutations: ", B, ".\n", sep="")
}
