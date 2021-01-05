
#' @title Simulations for p-Value Combinations
#' @description This function describes the behavior of a p-value combination in a given setting
#' through simulations.
#' @usage simulatePvals(type = "fisher", r = 1, rho = c(0,0.33,0.66,0.99), prop = c(0.1,0.2,0.3,0.4), m = 1000,
#'               B = 200, n = 50, alpha = 0.05, power = 0.8, rand = FALSE, truncFrom = alpha, truncTo = 1,
#'               nMax = 10000, nSim = 1000, silent = FALSE)
#' @param type transformation of p-values (\code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{vovk.wang})
#' @param r parameter for Vovk and Wang's p-value transformation.
#' @param rho vector of levels of equicorrelation between pairs of variables.
#' @param prop vector of proportions of non-null hypotheses.
#' @param m total number of variables.
#' @param B number of permutations, including the identity.
#' @param n number of observations.
#' @param alpha significance level.
#' @param power power of the one-sample t test.
#' @param rand logical, \code{TRUE} to compute p-values by permutation distribution.
#' @param truncFrom truncation parameter: values greater than \code{truncFrom} are truncated.
#' If \code{NULL}, p-values are not truncated.
#' @param truncTo truncation parameter: truncated values are set to \code{truncTo}.
#' If \code{NULL}, p-values are not truncated.
#' @param nMax maximum number of iterations.
#' @param nSim number of simulations per scenario.
#' @param silent logical, \code{FALSE} to print the summary.
#' @details For each scenario, the function generates
#' \code{nSim} matrices of simulated p-values through \code{\link{simulateMat}},
#' with seed equal to 1,...,\code{nSim}.
#' Then it determines a lower confidence bound for the number of true discoveries
#' within the set of non-null hypotheses.
#' @details The significance level \code{alpha} should be in the interval [1/\code{B}, 1).
#' @details Truncation parameters should be such that \code{truncTo} is not smaller than \code{truncFrom}.
#' @details A p-value \code{p} is transformed as following.
#' \itemize{
#' \item Edgington: \code{-p}
#' \item Fisher: \code{-log(p)}
#' \item Pearson: \code{log(1-p)}
#' \item Liptak: \code{-qnorm(p)}
#' \item Cauchy: \code{tan(0.5 - p)/p}
#' \item Vovk and Wang: \code{- sign(r)p^r}
#' }
#' @details Pearson's and Liptak's transformations produce infinite values in \code{1}.
#' For such transformations, \code{truncTo} is coerced to be not greater than \code{1 -  .Machine$double.eps}.
#' @return \code{simulatePvals} returns two lists, \code{results} and \code{summary}.
#' The list \code{results} contains, for each value of \code{prop}, the matrices
#' \code{TD}, \code{iterations} and \code{convergence}, which describe the results
#' for each value of \code{rho} and seed.
#' The list summary contains the matrices
#' \itemize{
#' \item \code{TDP}: matrix of the mean true discovery proportion for scenario
#' \item \code{iterations}: matrix of the mean number of iterations for scenario
#' \item \code{convergence}: logical, \code{TRUE} if the algorithm converges in all the simulations
#' }
#' @author Anna Vesely.
#' @examples
#' # 16 scenarions, 10 simulations each
#' sim <- simulatePvals(type="vovk.wang", r=-1, rho=c(0,0.33,0.66,0.99), prop = c(0.1,0.2,0.3,0.4), nSim=10)
#' @export


simulatePvals <- function(type="fisher", r=1, rho=c(0,0.33,0.66,0.99), prop = c(0.1,0.2,0.3,0.4), m=1000, B=200,
                          n=50, alpha=0.05, power=0.8, rand=FALSE, truncFrom=alpha, truncTo=1, nMax=10000,
                          nSim=1000, silent=FALSE){
  
  R <- length(rho)
  P <- length(prop)
  
  summary.TDP <- matrix(nrow=R, ncol=P)
  rownames(summary.TDP) <- rho
  colnames(summary.TDP) <- prop
  summary.it <- summary.TDP
  summary.conv <- summary.TDP
  
  
  results <- list()
  
  
  for(j in seq(P)){
    m1 <- round(prop[j]*m)
    S <- seq(m1)
    
    TD <- matrix(nrow=nSim, ncol=R)
    rownames(TD) <- seq(nSim)
    colnames(TD) <- rho
    it <- TD
    conv <- TD
    
    for(i in seq(R)){
      
      for(seed in seq(nSim)){
        G <- simulateMat(prop[j], m, B, rho[i], n, alpha, power=0.8, pvalues=TRUE, rand=FALSE, seed)
        res <- sumSome.pvals(G, S, alpha, truncFrom, truncTo, type, r, nMax)
        TD[seed,i] <- res$summary["TD"]
        it[seed,i] <- res$iterations
        conv[seed,i] <- (res$summary["TD"] == res$summary["maxTD"])*1
      }
      
      summary.TDP[i,j] <- round(mean(TD[,i])*100/m1, 1)
      summary.it[i,j] <- round(mean(it[,i]))
      summary.conv[i,j] <- all(conv[,i] == 1)*1
    }
    
    e <- list("prop"=prop[j], "TD"=TD, "iterations"=it, "convergence"=conv)
    results <- append(results, list(e))
  }
  
  summary <- list("TDP"=summary.TDP, "iterations"=summary.it, "convergence"=all(summary.conv==1))
  
  if(!silent){print(summary)}
  out <- list("results"=results, "summary"=summary)
  return(out)
}