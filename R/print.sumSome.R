


sumSome <- function(size, TD, maxTD, TDP, maxTDP, iterations){
  val <- list("size"=size, "TD"=TD, "maxTD"=maxTD, "TDP"=TDP, "maxTDP"=maxTDP, "iterations"=iterations)
  attr(val, "class") <- "sumSome"
  return(val)
}







sumBrain <- function(G, mask, alpha, truncFrom, truncTo){
  val <- list("statistics"=G, "mask"=mask, "alpha"=alpha, "truncFrom"=truncFrom, "truncTo"=truncTo)
  attr(val, "class") <- "sumBrain"
  return(val)
}







#' @export
print.sumSome <- function(obj){
  s <- as.character(obj$size)
  TD <- as.character(obj$TD)
  conv <- ifelse(obj$TD == obj$maxTD, "TRUE", "FALSE")
  
  cat("A sumSome object for the test of ", s, " hypotheses.\n", sep="")
  cat("True discoveries: ", TD, ".\n", sep="")
  cat("Convergence: ", conv, ".\n", sep="")
}



#' @export
print.sumBrain <- function(obj){
  f <- as.character(ncol(obj$statistics))
  cat("A sumBrain object containing brain imaging data with ", f, " voxels.\n", sep="")
}



