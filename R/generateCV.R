#' @title Vector of Critical Values for p-Value Combinations (Parametric)
#' @description Internal function.
#' It determines a vector of critical values for summed transformed p-values considering \code{1:m} hypotheses.
#' @usage generateCV(m, alpha, type, r, independence)
#' @param m number of hypotheses.
#' @param alpha significance level.
#' @param type p-value combination among \code{fisher}, \code{pearson}, \code{liptak}, \code{cauchy},
#' \code{harmonic}, \code{vovk.wang} (see details).
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param independence logical, \code{TRUE} to assume independence, \code{FALSE} for general dependence structure.
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
#' @details For \code{vovk.wang}, \code{r=-Inf} and \code{r=Inf} correspond to taking the minimum and the maximum
#' of the p-values, respectively.
#' @details Under general dependence, the test is defined only for \code{fisher}, \code{harmonic} and \code{vovk.wang}.
#' The latter always assumes general dependence.
#' @return \code{generateCV} returns a numeric vector of critical values of length \code{m}.
#' @author Xu Chen.
#' @noRd
#' @importFrom FMStable qEstable
#' @importFrom stats uniroot qchisq qcauchy


generateCV <- function(m, alpha, type, r, independence){
  
  
  # Fisher (independence)
  if(type=="fisher") {
    out <- qchisq(alpha, df=2*(1:m), lower.tail=FALSE)
    return(out)
  }
  
  
  # Pearson (independence)
  if(type=="pearson") {
    out <- -qchisq(alpha, df=2*(1:m))
    return(out)
  }
  
  
  # Liptak (independence)
  if(type=="liptak") {
    out <- sqrt(1:m)*qnorm(alpha, lower.tail=FALSE)
    return(out)
  }
  
  
  # Cauchy (independence)
  if(type=="cauchy"){
    out <- (1:m)*qcauchy(alpha, location=0, scale=1, lower.tail=FALSE)
    return(out)
  }
  
  
  
  # Harmonic mean (independence)
  if(type == "harmonic.ind"){
    out <- rep(1/alpha, m)
    if (m >= 2) {
      for (i in 2:m) {
        out[i] <- i * FMStable::qEstable(alpha, FMStable::setParam(alpha=1, location=(log(i)+1+digamma(1)-log(2/pi)), logscale=log(pi/2), pm=0), lower.tail=FALSE)
      }
    }
    return(out)
  }
  
  
  
  # GENERAL DEPENDENCE
  
  if(type=="vovk.wang"){
    
    # compute an
    if(m == 1){
      an <- 1
    }else{
      if(r == Inf){ # precise
        an <- rep(1,m)
      }else if(r == -Inf){ # precise
        an <- 1:m
      }else if(r >= m-1){ # precise
        an <- (1:m)^(1/r)
      }else if(r >= 1/(m-1) && r <= m-1){ # precise
        an <- c(1, rep((r+1)^(1/r),m-1))
      }else if(r == 0){ # precise
        an <- c(1, rep(exp(1),m-1))
        for (i in 2:min(16,m)){
          root  <- uniroot(function(x) log(1/x-(i-1)) - i + i^2*x, lower=0, upper=1/i-1e-8, tol=1e-30)$root
          an[i] <- 1/root * exp(-(i-1)*(1-i*root))
        }
      }else if(r == -1){ # precise
        an <- 1:m
        if (m >= 3) {
          for (i in (3:m)){
            root  <- uniroot(function(x) x^2 - i*(x+1)*log(x+1) + i*x, lower=i/2-1, upper=1e8, tol=1e-30)$root
            an[i] <- (root+i)^2/(root+1)/i
          }
        }
      }else if(r > -1){ # asymptotically precise
        an <- c(1, rep((r+1)^(1/r),m-1))
      }else if(r < -1){ # asymptotically precise
        an <- c(1, r/(r+1)*(2:m)^(1+1/r))
      }
    }
    
    # compute critical values
    if (r == 0) {
      out <- -(1:m) * log(alpha/an)
    } else {
      out <- -sign(r) * (1:m) * (alpha/an)^r
    }
    return(out)
  }
}

