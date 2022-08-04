#' @title Vector of Critical Values for p-Value Combinations (Parametric)
#' @description Internal function.
#' It determines a vector of critical values for summed transformed p-values considering \code{1:m} hypotheses.
#' @usage generateCV(m, type, alpha, r)
#' @param m number of hypotheses.
#' @param type p-value combination among \code{fisher}, \code{fisher.dep}, \code{pearson}, \code{liptak}, \code{cauchy}, \code{vovk.wang}, \code{harmonic.dep}, \code{harmonic.ind}
#' (see details).
#' @param alpha significance level.
#' @param r parameter for Vovk and Wang's p-value combination for general dependence.
#' @details The \code{type} determines the vector of critical values as following.
#' \itemize{
#' \item Fisher (independence, \code{fisher}): valid under independence, anti-conservative otherwise (Fisher, 1925)
#' \item Fisher (dependence, \code{fisher.dep}): valid under general dependence (Vovk and Wang, 2020)
#' \item Pearson (\code{pearson}): valid under independence, anti-conservative otherwise (Pearson, 1933)
#' \item Liptak (\code{liptak}): valid under independence, conservative or anti-conservative otherwise depending on the dependence structure among the tests (Liptak, 1958; unweighted version, same as Stouffer’s method (Stouffer et al., 1949))
#' \item Cauchy (\code{cauchy}): valid under independence and perfect dependence, approximately valid otherwise under bivariate normality assumption (Liu and Xie, 2020)
#' \item P-value combination for general dependence (\code{vovk.wang}): valid under general dependence (Vovk and Wang, 2020)
#' \item Harmonic mean (dependence, \code{harmonic.dep}): valid under general dependence (Vovk and Wang, 2020)
#' \item Harmonic mean (independence, \code{harmonic.ind}): valid under independence, anti-conservative otherwise (Wilson, 2019)
#' }
#' @return \code{generateCV} returns a numeric vector of critical values of length \code{m}.
#' @author Xu Chen, Anna Vesely.
#' @noRd
#' @importFrom FMStable qEstable
#' @importFrom stats uniroot qchisq qcauchy


generateCV <- function(m, type, alpha, r){

  type = match.arg(tolower(type), c("fisher",
                                    "fisher.dep",
                                    "pearson", 
                                    "liptak",
                                    "cauchy", 
                                    "vovk.wang",
                                    "harmonic.dep", 
                                    "harmonic.ind"))
  
  # (1) Fisher (independence)
  if(type=="fisher") {
    out <- qchisq(alpha, df=2*(1:m), lower.tail=FALSE)
    return(out)
  }
  
  # (2) Fisher (dependence)
  if(type=="fisher.dep") {
    an <- c(1, rep(exp(1),m-1))
    if (m >= 2) {
      for (i in 2:min(16,m)){
        root  <- uniroot(function(x) log(1/x-(i-1)) - i + i^2*x, lower=0, upper=1/i-1e-8, tol=1e-30)$root
        an[i] <- 1/root * exp(-(i-1)*(1-i*root))
      }
    }
    out <- -2 * (1:m) * log(alpha/an)
    return(out)
  }
  
  # (3) Pearson
  if(type=="pearson") {
    out <- -qchisq(alpha, df=2*(1:m))
    return(out)
  }
  
  # (4) Liptak
  if(type=="liptak") {
    out <- sqrt(1:m)*qnorm(alpha, lower.tail=FALSE)
    return(out)
  }
  
  # (5) Cauchy
  if(type=="cauchy"){
    out <- (1:m)*qcauchy(alpha, location=0, scale=1, lower.tail=FALSE)
    return(out)
  }
  
  # (6) P-value combination for general dependence
  if(type=="vovk.wang"){
    
    # compute an
    if (m == 1) {
      an <- 1
    } else if (m >= 2) {
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
  
  # (7) Harmonic mean (dependence)
  if(type == "harmonic.dep"){
    an <- 1:m
    if (m >= 3) {
      for (i in 3:m) {
        root  <- uniroot(function(x) x^2-i*(x+1)*log(x+1)+i*x, lower=i/2-1, upper=1e8, tol=1e-30)$root
        an[i] <- (root+i)^2/(root+1)/i
      }
    }
    out <- (1:m)*an/alpha
    return(out)
  }
  
  # (8) Harmonic mean (independence)
  if(type == "harmonic.ind"){
    out <- rep(1/alpha, m)
    if (m >= 2) {
      for (i in 2:m) {
        out[i] <- i * FMStable::qEstable(alpha, FMStable::setParam(alpha=1, location=(log(i)+1+digamma(1)-log(2/pi)), logscale=log(pi/2), pm=0), lower.tail=FALSE)
        #out[i] <- i * FMStable::qEstable(1-alpha, FMStable::setParam(alpha=1, location=(log(i)+1+digamma(1)-log(2/pi)), logscale=log(pi/2), pm=0))
      }
    }
    return(out)
  }
  
}