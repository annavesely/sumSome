#' @title True Discovery Guarantee for p-Value Combinations - Parametric
#' @description This function uses p-values to determine confidence bounds
#' for the number of true discoveries, the true discovery proportion
#' and the false discovery proportion within a set of interest.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage sumPvalsPar(g, S = NULL, alpha = 0.05, type = "vovk.wang", r = 0, independence = NULL)
#' @param g numeric vector of p-values.
#' @param S vector of indices for the variables of interest (if not specified, all variables).
#' @param alpha significance level.
#' @param type p-value combination among \code{fisher}, \code{pearson}, \code{liptak}, \code{cauchy},
#' \code{harmonic}, \code{vovk.wang} (see details).
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param independence logical, \code{TRUE} to assume independence, \code{FALSE} for general dependence structure.
#' If not specified, it is set to \code{FALSE} for \code{vovk.wang}, and \code{TRUE} otherwise.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Fisher: \code{-2log(p)} (Fisher, 1925)
#' \item Pearson: \code{2log(1-p)} (Pearson, 1933)
#' \item Liptak: \code{qnorm(1-p)} (Liptak, 1958; Stouffer et al., 1949)
#' \item Cauchy: \code{tan[(0.5-p)pi]} with \code{pi=3.142} (Liu and Xie, 2020)
#' \item Harmonic mean: \code{1/p} (Wilson, 2019)
#' \item Vovk and Wang: \code{p^r} (\code{log(p)} for \code{r}=0) (Vovk and Wang, 2020)
#' }
#' An error message is returned if the transformation produces infinite values.
#' @details For Vovk and Wang, \code{r=-Inf} corresponds to the minimum p-value, \code{r=Inf} to the maximum p-value,
#' \code{r=0} to Fisher, and \code{r=-1} to the harmonic mean.
#' @details Under independence, for Vovk and Wang the test is defined only
#' for \code{r=0} and \code{r=1}. Under general dependence, the test is defined only for
#' Fisher, the harmonic mean and Vovk and Wang.
#' @details For combinations that are not implemented, if the vector of critical values is known
#' the method can be applied through \code{\link{sumStatsPar}}.
#' Please contact us to implement other known vectors of critical values that do not currently appear.
#' @return \code{sumPvalsPar} returns an object of class \code{sumObj}, containing
#' \itemize{
#' \item \code{total}: total number of variables (length of \code{g})
#' \item \code{size}: size of \code{S}
#' \item \code{alpha}: significance level
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
#' \item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
#' \item \code{iterations}: number of iterations of the algorithm (\code{NULL})
#' }
#' @author Xu Chen.
#' @examples
#' # generate vector of p-values for 5 variables
#' g <- as.vector(simData(prop = 0.6, m = 5, B = 1, alpha = 0.4, seed = 42))
#' 
#' # subset of interest (variables 1 and 2)
#' S <- c(1,2)
#'  
#' # create object of class sumObj
#' # combination: harmonic mean under general dependence
#' res <- sumPvalsPar(g, S, alpha = 0.4, type = "harmonic", independence = FALSE)
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
#' @references
#' Goeman J. J. and Solari A. (2011). Multiple testing for exploratory research. Statistical Science, doi: 10.1214/11-STS356.
#' 
#' Tian J., Chen X., Katsevich E., Goeman J. J. and Ramdas A. (2022). Large-scale simultaneous inference under dependence. Scandinavian Journal of Statistics, doi: 10.1111/sjos.12614.
#' @seealso
#' True discovery guarantee using generic statistics (parametric): \code{\link{sumStatsPar}}
#' 
#' Access a \code{sumObj} object: \code{\link{discoveries}}, \code{\link{tdp}}, \code{\link{fdp}}
#' @export

sumPvalsPar <- function(g, S=NULL, alpha=0.05, type="vovk.wang", r=0, independence=NULL){
  
  # change lines 95 and 104
  
  if(!is.vector(g) || !is.numeric(g) || !all(is.finite(g))){stop("g must be a vector of finite numbers")}
  if(length(g)==0){stop("g must be a vector of finite numbers")}
  if(!all(g >= 0) || !all(g <= 1)){stop("g must be a vector of p-values")}
  
  if(is.null(S)){S <- seq(length(g))}
  if(!is.vector(S) || !is.numeric(S) || !all(is.finite(S))){stop("S must be a vector of finite integers")}
  if(!all(floor(S)==S)){stop("S must be a vector of finite integers")}
  if(!all(S >= 0) || !all(S <= length(g))){stop("S must contain indices between 1 and the total number of variables")}
  S <- unique(S)
  
  if(!is.numeric(alpha) || !is.finite(alpha)){stop("alpha must be a number in (0,1)")}
  if(alpha <= 0 || alpha >= 1){stop("alpha must be a number in (0,1)")}
  
  type = match.arg(tolower(type), c("fisher", "pearson", "liptak", "cauchy", "harmonic", "vovk.wang"))
  
  if(!is.numeric(r)){stop("r must be a real number")}
  
  # if independence is not specified, we assume it is TRUE unless the type is Vovk and Wang
  if(is.null(independence)){
    independence <- (type != "vovk.wang")
  }
  
  # switch to Fisher & Harmonic mean for independence
  if(type == "vovk.wang" && independence){
    if(r == 0){type <- "fisher"}
    else if(r == -1){type <- "harmonic"}
    else{
      stop("The method is not implemented under independence for this p-value combination.
  Use independence = FALSE to assume general dependence.
  Otherwise, refer to sumStatsPar() to manually input a vector of critical values.")
    }
  }
  
  # switch to Vovk and Wang for general dependence
  if(type != "vovk.wang" && !independence){
    if(type == "fisher"){r <- 0}
    else if(type == "harmonic"){r <- -1}
    else{
      stop("The method is not implemented under general dependence for this p-value combination.
  Use independence = TRUE to assume independence.
  Otherwise, refer to sumStatsPar() to manually input a vector of critical values.")
    }
    type <- "vovk.wang"
  }
  
  # transform p-values
  if(type=="fisher"){g <- -2*log(g)}
  else if(type=="pearson"){g <- 2*log(1-g)}
  else if(type=="liptak"){g <- qnorm(1-g)}
  else if(type=="cauchy") {g <- tan((0.5-g)*pi)}
  else if (type=="harmonic"){g <- 1/g}
  else if(type=="vovk.wang"){
    if(r==0){g <- -log(g)}
    else if(r==Inf){g <- max(g)}
    else if(r==-Inf){g <- min(g)}
    else{g <- -sign(r) * g^r}
  }
  
  if(!all(is.finite(g))){stop("Transformation produced infinite values")}
  
  # generate critical vector
  cvs <- generateCV(length(g), alpha, type, r, independence)
  
  out <- sumTestPar(g, S, alpha, cvs)
  return(out)
}