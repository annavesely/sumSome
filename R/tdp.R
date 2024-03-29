#' @title Confidence Bound for the TDP
#' @description This function determines a lower confidence bound for the true discovery proportion
#' within a set of interest. The bound remains valid under post-hoc selection.
#' @usage tdp(object)
#' @param object an object of class \code{sumObj}, as returned by
#' the functions \code{\link{sumStats}} and \code{\link{sumPvals}}.
#' @return \code{tdp} returns a lower (1-\code{alpha})-confidence bound
#' for the true discovery proportion in the set.
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
#' 
#' Lower confidence bound for the number of true discoveries: \code{\link{discoveries}}
#' 
#' Upper confidence bound for the FDP: \code{\link{fdp}}
#' @export

tdp <- function(object){
  UseMethod("tdp")
}



#' @rdname tdp
#' @export

tdp.sumObj = function(object) {
  return(object$TD/object$size)
}