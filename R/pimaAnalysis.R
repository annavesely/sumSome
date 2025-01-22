#' True discovery guarantee in multiverse analysis
#' @description This function uses permutation score statistics/p-values to determine a true discovery guarantee for multiverse analysis,
#' when studying one or more parameters of interest within a multiverse of models.
#' It computes confidence bounds for the number of true discoveries and the true discovery proportion overall or within different groups.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage pimaAnalysis(pima, by = NULL, type = "sum", r = 0, alpha = 0.05, ...) 
#' @param pima an object of class \code{pima}, as obtained from the function \code{\link[pima]{pima}}.
#' @param by name of grouping element among \code{Coeff} and \code{Model}. If not specified, all coefficients of interest
#' in all models are considered together.
#' @param type combining function: \code{sum} uses the sum of scores as in \code{\link[sumSome]{sumStats}},
#' while different p-value combinations are defined as in \code{\link[sumSome]{sumPvals}}.
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param alpha significance level.
#' @param ... further parameters of \code{sumStats} or \code{sumPvals} (truncation parameters and maximum number of iterations of the algorithm).
#' @details In the default \code{by = NULL}, the procedure computes lower confidence bounds for the
#' number/proportion of significant effects (non-null coefficients) among all.
#' Other inputs of the argument \code{by} return analogous bounds, defined by coefficient (\code{"Coeff"}) or by model (\code{"Model"}).
#' @details If truncation parameters are not specified, scores/p-values are not truncated.
#' @return Returns a data frame containing a summary for each subset:
#' \itemize{
#' \item \code{size}: number of considered coefficients
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the true discovery proportion
#' }
#' @references
#' Girardi, Vesely, Lakens, Altoè, Pastore, Calcagnì, Finos (2024). Post-selection Inference in Multiverse Analysis (PIMA): An Inferential Framework Based on the Sign Flipping Score Test. Psychometrika, doi: 10.1007/s11336-024-09973-6.
#' 
#' Vesely, Finos, Goeman (2023). Permutation-based true discovery guarantee by sum tests. Journal of the Royal Statistical Society, Series B (Statistical Methodology), doi: 10.1093/jrsssb/qkad019
#' @examples
#' # COPIARE PARTE DELL'ESEMPIO DELLA FUNZIONE pima()
#' @export

pimaAnalysis <- function(pima, by = NULL, type = "sum", r = 0, alpha = 0.05, ...){
  
  type <- match.arg(tolower(type), c("sum", "fisher", "pearson", "liptak",  "edgington", "cauchy", "harmonic", "vovk.wang"))
  
  extra_args <- list(...)
  truncFrom <- extra_args$truncFrom %||% NULL
  truncTo <- extra_args$truncTo %||% NULL
  
  if(is.null(by)){
    S = list(Overall = 1:ncol(obj$Tspace))
  }else{
    smr = apply(obj$summary_table[, by, drop = FALSE], 1, paste, collapse = ".")
    uniq_nm = unique(smr)
    S = lapply(uniq_nm, function(nm) which(smr == nm))
    names(S) = uniq_nm
  }
  
  G <- abs(as.matrix(res$Tspace))
  alpha <- 1 - conf_level
  
  if(type == "sum"){
    
    get_disc <- function(x){
      tmp <- do.call(sumSome::sumStats, c(
        list(G=G, S=S[[x]], alternative="greater", alpha=alpha), 
        extra_args
      ))
      return(discoveries(tmp))
    }
  }else{
    
    P <- flip::t2p(G, obs.only = FALSE)
    
    get_disc <- function(x){
      tmp <- do.call(sumSome::sumPvals, c(
        list(G=G, S=S[[x]], alpha=alpha, type = type, r = r), 
        extra_args
      ))
      return(discoveries(tmp))
    }
  }
  
  out <- data.frame(
    size = sapply(1:length(S), FUN = function(i) length(S[[i]])),
    TD = sapply(1:length(S), FUN=get_disc)
  )
  
  out$TDP <- round(out$TD/out$size, 3)
  names(out) <- names(S)
  
  return(out)
}

