#' @title Vector of Critical Values for p-Value Combinations (Parametric)
#' @description Internal function.
#' It determines a vector of critical values for summed transformed p-values considering \code{1:m} hypotheses.
#' @usage generateCV(m, type, alpha)
#' @param m number of hypotheses.
#' @param type p-value combination among \code{harmonic.dep}, \code{harmonic.ind}, \code{fisher.ind}, \code{cauchy.ind}
#' (see details).
#' @param alpha significance level.
#' @details The \code{type} determines the vector of critical values as following.
#' \itemize{
#' \item Harmonic mean under dependence: valid under general dependence (Vovk and Wang, 2020)
#' \item Harmonic mean under independence: valid under independence, anti-conservative under general dependence (Wilson, 2019)
#' \item Fisher under independence:  valid under independence, anti-conservative under general dependence (Fisher, 1925)
#' \item Cauchy: valid under independence and perfect dependence, approximately valid under general dependence (Liu and Xie, 2020)
#' }
#' @return \code{generateCV} returns a numeric vector of critical values of length \code{m}.
#' @author Xu Chen, Jelle Goeman.
#' @noRd


generateCV <- function(m, type, alpha){
  
  type = match.arg(tolower(type), c("harmonic.dep", "harmonic.ind", "fisher.ind", "cauchy"))
  
  # (1) HMP - Vovk & Wang's (valid for general dep.)
  if(type == "harmonic.dep"){
    an <- numeric(m)+1
    an[2] <- 2
    for (i in 3:m) {
      root  <- uniroot(function(x) x^2-i*(x+1)*log(x+1)+i*x, lower=i/2-1, upper=1e8, tol=1e-9)$root
      an[i] <- (root+i)^2/(root+1)/i
    }
    out <- (1:m)*an/alpha
    return(out)
  }
  
  # (2) HMP - Wilson's (valid for ind., but anti-conservative for dep.)
  if(type == "harmonic.ind"){
    out <- numeric(m)+alpha
    for (i in 2:m) {
      out[i] <- i*FMStable::qEstable(1-alpha,setParam(alpha=1,location=(log(i)+1+digamma(1)-log(2/pi)),logscale=log(pi/2),pm=0))
    }
    return(out)
  }
  
  # (3) Fisher combination test (valid for ind., but anti-conservative for dep.)
  if(type=="fisher.ind") {
    out <- qchisq(1-alpha, df=2*(1:m))
    return(out)
  }
  
  # (4) Cauchy combination test (valid for ind. & perfect dep., and approx. valid for dep.)
  if(type=="cauchy"){
    out <- (1:m)*qcauchy(alpha, location=0, scale=1, lower.tail=FALSE)
    return(out)
  }
}