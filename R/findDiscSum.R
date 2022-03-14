#' @title Find Discoveries Using Summed Statistics
#' @name findDiscSum
#' @aliases findDiscSum
#' @description Based on a newly developed linear-time shortcut, this function determines lower confidence bounds for the number of true discoveries for a set of interest.
#' @usage findDiscSum(stats, cvs, ix = NULL, alpha = 0.05)
#' @param stats Numeric vector of \code{m} statistic values. If \code{cvs} is numeric, \code{stats} should be transformed statistics; if \code{cvs} is of type "character", \code{stats} should contain original p-values.
#' @param cvs Critical values for sum tests. Can be either a (numeric) \code{m}-vector of critical values for summed transformed statistics considering \code{1:m} hypotheses, or a (character) vector specifying the type of p-value combination test ("hmp.dep", "hmp.ind", "fisher.ind" or "cauchy").
#' @param ix Indices of subset of interest \code{S}. If \code{ix = NULL}, it is assumed that the full set is considered.
#' @param alpha Significance level.
#' @return \code{findDiscSum} returns a lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}.
#' @author Xu Chen, Jelle Goeman.
#' @examples
#' # generate matrix of p-values for 5 variables and 10 permutations
#' G <- simData(prop = 0.6, m = 5, B = 10, alpha = 0.4, p = TRUE, seed = 42)
#' 
#' # compute lower confidence bound for the number of true discoveries for the whole set using harmonic mean p-value
#' res <- findDiscSum(stats = G[1,], cvs = "hmp.ind")
#' 
#' @import FMStable
#' @export

findDiscSum <- function(stats, cvs, ix=NULL, alpha=0.05) {
  # check ix
  if (is.null(ix)) {ix <- seq(m)}
  m <- length(stats)
  s <- length(ix)
  
  if (m>0) {
    # check cvs (input critical values)
    if (missing(cvs)) {
      stop("Please specify critical values.")
    } else if (is.character(cvs)) {
      cvs <- match.arg(tolower(cvs), c("hmp.dep","hmp.ind","fisher.ind","cauchy"))
      if (cvs=="hmp.dep" || cvs=="hmp.ind") {
        stats <- 1/stats
      } else if (cvs=="fisher.ind") {
        stats <- -2*log(stats)
      } else if (cvs=="cauchy") {
        stats <- tan((0.5-stats)*pi)
      }
      cvs <- generateCV(m, type=cvs, alpha=alpha)
    } else if (!is.numeric(cvs)) {
      stop("Critical values should be numeric.")
    }
    
    notix <- (1:m)[-ix]
    u0 <- max(stats)
    v0 <- min(stats)
    u <- c(v0, sort(stats[ix]), rep(u0+1, m-s+1))
    v <- c(v0, sort(stats[notix]), rep(u0+1, s+1))
    
    out <- findDiscov_sum(s, m, u, v, cvs)
  } else {
    out <- 0
  }
  
  return(out)
}

# Compute critical values for summed transformed p-values
generateCV <- function(m, type="hmp.dep", alpha=0.05) {

  if (type=="hmp.dep") {
    # (1) HMP - Vovk & Wang's (valid for general dep.)
    an <- numeric(m)+1
    an[2] <- 2
    for (i in 3:m) {
      root  <- uniroot(function(x) x^2-i*(x+1)*log(x+1)+i*x, lower=i/2-1, upper=1e8, tol=1e-9)$root
      an[i] <- (root+i)^2/(root+1)/i
    }
    out <- (1:m)*an/alpha
  } else if (type=="hmp.ind") {
    # (2) HMP - Wilson's (valid for ind., but anti-conservative for dep.)
    out <- numeric(m)+alpha
    for (i in 2:m) {
      out[i] <- i*FMStable::qEstable(1-alpha,setParam(alpha=1,location=(log(i)+1+digamma(1)-log(2/pi)),logscale=log(pi/2),pm=0))
    }
  } else if (type=="fisher.ind") {
    # (3) Fishier combination test (valid for ind., but anti-conservative for dep.)
    out <- qchisq(1-alpha, df=2*(1:m))
  } else if (type=="cauchy") {
    # (4) Cauchy combination test (valid for ind. & perfect dep., and approx. valid for dep.)
    out <- (1:m)*qcauchy(alpha, location=0, scale=1, lower.tail=FALSE)
  } else {
    stop("Please select an available sum test type.")
  }
  
  return(out)
}

