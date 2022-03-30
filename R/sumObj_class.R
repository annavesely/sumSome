#' @title sumObj Class
#' @description Internal function.
#' It defines a \code{sumObj} object, storing information on the closed testing procedure on a set of interest.
#' @usage sumObj(total, size, alpha, TD, maxTD, iterations)
#' @param total total number of variables.
#' @param size size of the set of interest.
#' @param alpha significance level.
#' @param TD lower (1-\code{alpha})-confidence bound for the number of true discoveries in the set.
#' @param maxTD maximum value of \code{TD} that could be found under convergence of the algorithm.
#' @param iterations number of iterations of the algorithm.
#' @return \code{sumObj} returns an object of class \code{sumObj}, containing
#' \code{total}, \code{size}, \code{alpha}, \code{TD}, \code{maxTD} and \code{iterations}.
#' @author Anna Vesely.
#' @noRd

sumObj <- function(total, size, alpha, TD, maxTD, iterations){
  val <- list("total"=total, "size"=size, "alpha"=alpha, "TD"=TD, "maxTD"=maxTD, "iterations"=iterations)
  attr(val, "class") <- "sumObj"
  return(val)
}





#' @export

print.sumObj <- function(x, ...){
  s <- ifelse(x$size==1, "1 hypothesis", paste(as.character(x$size), " hypotheses", sep=""))
  f <- as.character(x$total)
  alpha <- as.character(x$alpha)
  
  cat("A sumObj object for closed testing on ", s, " out of ", f, ", with significance level ", alpha, ".\n", sep="")
  cat("Use discoveries(), tdp() or fdp() to access this object.\n", sep="")
}





#' @export

summary.sumObj <- function(object, ...){
  TD <- ifelse(object$TD==1, "1 discovery", paste(as.character(object$TD), " discoveries", sep=""))
  TDP <- as.character(round(object$TD/object$size, 2))
  conf <- as.character(1 - object$alpha)
  conv <- ifelse(object$TD == object$maxTD, "converged", "did not converge")
  
  print(object)
  cat("\n", sep="")
  cat("With ", conf, " confidence: at least ", TD, " (a proportion of ", TDP, ").\n", sep="")
  
  if(!is.null(object$iterations)){
    iter <- ifelse(object$iterations==1, "1 iteration", paste(as.character(object$iterations), " iterations", sep=""))
    cat("Algorithm ", conv, " after ", iter, ".\n", sep="")
  }
}










