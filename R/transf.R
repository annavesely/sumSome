#' @title Statistics Transformation 
#' @description Internal function. It truncates and transforms a matrix of statistics.
#' @usage transf(G, truncFrom, truncTo, option, r)
#' @param G matrix of statistics
#' @param truncFrom truncation parameter: values less extreme than \code{truncFrom} are truncated
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}
#' @param option direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}),
#' or transformation (\code{squares} for generic statistics,
#' and \code{edgington}, \code{fisher}, \code{pearson}, \code{liptak}, \code{cauchy}, \code{vovk.wang} for p-values)
#' @param r parameter for the Vovk and Wang p-value transformation
#' @details The transformation is determined by the \code{option} as following. A generic value \code{x} is transformed to
#' \code{x} for \code{greater}, \code{-x} for \code{lower}, \code{|x|} for \code{two.sided},
#' \code{x^2} for \code{squares}, \code{-x} for \code{edgington}, \code{-log(x)} for \code{fisher},
#' \code{log(1-x)} for \code{pearson}, \code{-qnorm(x)} for \code{liptak}, \code{tan(0.5 - x)/x} for \code{cauchy},
#' and \code{- sign(r)x^r} for \code{vovk.wang}.
#' @details Pearson's and Liptak's transformations produce infinite values in \code{1}.
#' Hence when the \code{option} is \code{pearson} or \code{liptak}, the truncation parameter \code{truncTo} is coerced to be
#' not greater than \code{1 -  .Machine$double.eps}.
#' @return \code{transf} returns a list containing the truncated and transformed matrix \code{G}, and the transformed truncation parameters \code{truncFrom} and \code{truncTo}.
#' Those are determined so that the most extreme values are always the greatest.
#' @author Anna Vesely
#' @export



transf <- function(G, truncFrom, truncTo, option, r){
  
  pvalues <- !(option %in% c("greater", "lower", "two.sided", "squares"))
  
  if(!is.matrix(G) || !is.numeric(G) || !is.finite(G)){stop("G must be a matrix of finite numbers")}
  if(pvalues && (!all(G >= 0) || !all(G <= 1))){stop("G must be a matrix of pvalues")}
  
  if(!is.numeric(r) || !is.finite(r)){stop("r must be a finite number")}
  if(option=="vovk.wang" && r==0){option <- "fisher"}
  
  if(option == "lower" || option == "edgington"){G <- - G}
  else if(option == "two.sided"){G <- abs(G)}
  else if(option == "squares"){G <- G^2}
  else if(option == "fisher"){G <- -log(G)}
  else if(option == "pearson"){G <- log(1-G)}
  else if(option == "liptak"){G <- -qnorm(G)}
  else if(option == "cauchy"){G <- tan(0.5-G)/G}
  else if(option == "vovk.wang"){G <- - sign(r) * G^r}
  
  truncation <- (!is.null(truncFrom) && !is.null(truncTo))
  
  if(truncation){
    if(!is.numeric(truncFrom) || !is.finite(truncFrom)){stop("truncFrom must be a finite number")}
    if(!is.numeric(truncTo) || !is.finite(truncTo)){stop("truncTo must be a finite number")}
    if(pvalues && (truncFrom < 0 || truncFrom > 1)){stop("truncFrom must be a number in [0,1]")}
    if(pvalues && (truncTo < 0 || truncTo > 1)){stop("truncTo must be a number in [0,1]")}
    
    eps <- .Machine$double.eps
    
    if(option == "lower" || option == "edgington"){
      truncFrom <- - truncFrom
      truncTo <- - truncTo
    }else if(option == "two.sided"){
      truncFrom <- abs(truncFrom)
      truncTo <- abs(truncTo)
    }else if(option == "squares"){
      truncFrom <- truncFrom^2
      truncTo <- truncTo^2
    }else if(option == "fisher"){
      truncFrom <- - log(truncFrom)
      truncTo <- - log(truncTo)
    }else if(option == "pearson"){
      truncFrom <- max(log(1-truncFrom), log(eps))
      truncTo <- max(log(1-truncTo), log(eps))
    }else if(option == "liptak"){
      truncFrom <- max(qnorm(1-truncFrom), qnorm(eps))
      truncTo <- max(qnorm(1-truncTo), qnorm(eps))
    }else if(option == "cauchy"){
      truncFrom <- tan(0.5-truncFrom)/truncFrom
      truncTo <- tan(0.5-truncTo)/truncTo
    }else if(option == "vovk.wang"){
      truncFrom <- - sign(r) * truncFrom^r
      truncTo <- - sign(r) * truncTo^r
    }
    
    if(truncTo > truncFrom){stop("Invalid truncation parameters: truncTo must be more extreme than truncFrom")}
    
    for(i in seq(ncol(G))){
      for(b in seq(nrow(G))){
        if(G[b,i] < truncFrom){G[b,i] <- truncTo}
      }
    }
  }
  
  if(!all(is.finite(G))){stop("Transformation produced infinite values")}
  out <- list("G"=G, "truncFrom"=truncFrom, "truncTo"=truncTo)
  return(out)
}