#' @title True discovery guarantee
#' @name true-discovery-guarantee
#' @rdname true-discovery-guarantee
#' @description These functions determine a lower confidence bound for the number of true discoveries,
#' a lower confidence bound for the true discovery proportion (TDP),
#' and an upper confidence bound for the false discovery proportion (FDP)
#' within a set of interest. The bounds remain valid under post-hoc selection.
#' @usage discoveries(object)
#' @param object an object of class \code{sumObj}, as returned by
#' the functions \code{\link{sumStats}} and \code{\link{sumPvals}}.
#' @return \code{discoveries}, \code{tdp} and \code{fdp} return a (1-\code{alpha})-confidence bound for the corresponding quantity in the subset.
#' @author Anna Vesely.
#' @examples
#' # generate matrix of p-values for 5 variables and 10 permutations
#' G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, seed = 42)
#' 
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#'  
#' # create object of class sumObj
#' # combination: harmonic mean (Vovk and Wang with r = -1)
#' res <- sumPvals(G, S, alpha = 0.4, r = -1)
#' res
#' summary(res)
#' 
#' # lower confidence bound for the number of true discoveries in S
#' discoveries(res)
#' 
#' # lower confidence bound for the true discovery proportion in S
#' tdp(res)
#' 
#' # upper confidence bound for the false discovery proportion in S
#' fdp(res)
#' @seealso
#' Create a \code{sumObj} object: \code{\link{sumStats}}, \code{\link{sumPvals}}
#' @export

discoveries <- function(object){
  UseMethod("discoveries")
}


#' @export

discoveries.sumObj = function(object) {
  return(object$TD)
}



#' @rdname true-discovery-guarantee
#' @usage tdp(object)
#' @export

tdp <- function(object){
  UseMethod("tdp")
}

#' @export

tdp.sumObj = function(object) {
  return(object$TD/object$size)
}


#' @rdname true-discovery-guarantee
#' @usage fdp(object)
#' @export

fdp <- function(object){
  UseMethod("fdp")
}


#' @export

fdp.sumObj = function(object) {
  return(1- object$TD/object$size)
}
