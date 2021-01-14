#' @title sumSome Class
#' @description Internal function.
#' It defines a \code{sumSome} object, storing information on the closed testing procedure on a set of interest.
#' @usage sumSome(total, size, alpha, TD, maxTD, iterations)
#' @param total total number of variables.
#' @param size size of the set of interest.
#' @param alpha significance level.
#' @param TD lower (1-\code{alpha})-confidence bound for the number of true discoveries in the set.
#' @param maxTD maximum value of \code{TD} that could be found under convergence of the algorithm.
#' @param iterations number of iterations of the algorithm.
#' @return \code{sumSome} returns an object of class \code{sumSome}, containing
#' \code{total}, \code{size}, \code{alpha}, \code{TD}, \code{maxTD} and \code{iterations}.
#' @author Anna Vesely.
#' @keywords Internal

sumSome <- function(total, size, alpha, TD, maxTD, iterations){
  val <- list("total"=total, "size"=size, "alpha"=alpha, "TD"=TD, "maxTD"=maxTD, "iterations"=iterations)
  attr(val, "class") <- "sumSome"
  return(val)
}





#' @export

print.sumSome <- function(sumSome){
  s <- ifelse(sumSome$size==1, "1 hypothesis", paste(as.character(sumSome$size), " hypotheses", sep=""))
  f <- as.character(sumSome$total)
  alpha <- as.character(sumSome$alpha)
  
  cat("A sumSome object for closed testing on ", s, " out of ", f, ", with significance level ", alpha, ".\n", sep="")
  cat("Use discoveries(), tdp() or fdp() to access this object.\n", sep="")
}





#' @export

summary.sumSome <- function(sumSome){
  TD <- ifelse(sumSome$TD==1, "1 discovery", paste(as.character(sumSome$TD), " discoveries", sep=""))
  TDP <- as.character(round(sumSome$TD/sumSome$size, 2))
  conf <- as.character(1 - sumSome$alpha)
  conv <- ifelse(sumSome$TD == sumSome$maxTD, "converged", "did not converge")
  iter <- as.character(sumSome$iterations)
  
  print(sumSome)
  cat("\n", sep="")
  cat("With ", conf, " confidence: at least ", TD, " (a proportion of ", TDP, ").\n", sep="")
  cat("Algorithm ", conv, " after ", iter, " iterations.\n", sep="")
}










