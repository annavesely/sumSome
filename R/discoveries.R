
#' A description of MyHappyFunction
#'
#' A details of MyHappyFunction
#'
#' @title MyHappyFunction: The my happy function
#' @param x numeric number
#' @param ... other arguments
#' @examples
#' a <- 1
#' class(a) <- "lm"
#' MyHappyFunction(a)
#' @export
discoveries <- function(obj){
  UseMethod("discoveries")
}

#' @export
discoveries.sumSome = function(obj) {
  out <- list("TD"=obj$TD, "maxTD"=obj$maxTD)
  return(out)
}
