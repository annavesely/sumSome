#' True discovery guarantee in multiverse analysis
#' @description This function uses permutation statistics/p-values to determine a true discovery guarantee for multiverse analysis,
#' when studying one or more parameters of interest within a multiverse of models.
#' It computes confidence bounds for the number of true discoveries and the true discovery proportion overall or within different groups.
#' The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
#' @usage pimaAnalysis(obj, by = NULL, type = "sum", r = 0, alpha = 0.05, ...) 
#' @param obj an object of class \code{jointest}, as obtained from the functions \code{pima} (package \strong{pima}) or \code{join_flipscores} (\strong{jointest}).
#' @param by name of grouping element among \code{Coeff} and \code{Model}. If not specified, all coefficients of interest
#' in all models are considered together.
#' @param type combining function: \code{sum} uses the sum of test statistics as in \code{sumStats},
#' while different p-value combinations are defined as in \code{sumPvals} (\code{edgington}, \code{fisher}, \code{pearson}, \code{liptak},
#' \code{cauchy}, \code{harmonic}, \code{vovk.wang}).
#' @param r parameter for Vovk and Wang's p-value combination.
#' @param alpha significance level.
#' @param ... further parameters of \code{sumStats} or \code{sumPvals} (truncation parameters and maximum number of iterations of the algorithm).
#' @details In the default \code{by = NULL}, the procedure computes lower confidence bounds for the
#' number/proportion of significant effects (non-null coefficients) among all.
#' Other inputs of the argument \code{by} return analogous bounds, defined by coefficient (\code{"Coeff"}) or by model (\code{"Model"}).
#' While the bounds are simultaneous over all possible groupings,
#' the combining function \code{type} should be fixed in advance.
#' @details If truncation parameters are not specified among the further parameters,
#' statistics/p-values are not truncated.
#' @details More generically, \code{obj} can be any list containing:
#' \itemize{
#' \item \code{Tspace}: data frame of statistics, where columns correspond to variables,
#' and rows to data transformations (e.g. permutations). The first transformation is the identity.
#' \item \code{summary_table}: summary data frame where rows correspond to variables.
#' }
#' In this framework, the grouping element \code{by} is the name of a column of \code{summary_table}.
#' @return \code{pimaAnalysis} returns a data frame containing a summary for each subset:
#' \itemize{
#' \item \code{size}: number of considered coefficients
#' \item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of significant effects
#' \item \code{TDP}: lower (1-\code{alpha})-confidence bound for the proportion of significant effects
#' }
#' @references
#' Girardi P., Vesely A., Lakens D., Altoè G., Pastore M., Calcagnì A., and Finos L. (2024). Post-selection Inference in Multiverse Analysis (PIMA): An Inferential Framework Based on the Sign Flipping Score Test. Psychometrika, doi: 10.1007/s11336-024-09973-6.
#' 
#' Vesely A., Finos L., and Goeman J. J. (2023). Permutation-based true discovery guarantee by sum tests. Journal of the Royal Statistical Society, Series B (Statistical Methodology), doi: 10.1093/jrsssb/qkad019.
#' @examples
#' # generate matrix of statistics for 2 coefficients X and Z within 3 models
#' G <- simData(prop = 0.6, m = 6, B = 50, alpha = 0.4, p = FALSE, seed = 42)
#' colnames(G) <- rep(c("X","Z"),3)
#'  
#' # summary table
#' summary_table <- data.frame(
#'   Model = rep(c("mod1","mod2","mod3"), each=2),
#'   Coeff = colnames(G)
#' )
#' 
#' # list of Tspace and summary_table
#' obj <- list(Tspace = as.data.frame(G), summary_table = summary_table)
#' 
#' # significant effects overall (sum of test statistics)
#' pimaAnalysis(obj, alpha = 0.4)
#' 
#' # significant effects by coefficient (sum of test statistics)
#' pimaAnalysis(obj, by = "Model", alpha = 0.4)
#' 
#' # significant effects by model (Fisher's combination of p-values)
#' pimaAnalysis(obj, by = "Coeff", type = "fisher", alpha = 0.4)
#' 
#' @seealso
#' True discovery guarantees: \code{\link{sumStats}}, \code{\link{sumPvals}}
#' @export

pimaAnalysis <- function(obj, by = NULL, type = "sum", r = 0, alpha = 0.05, ...){
  
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
  
  G <- abs(as.matrix(obj$Tspace))
  
  if(type == "sum"){
    
    get_disc <- function(x){
      tmp <- do.call(sumSome::sumStats, c(
        list(G=G, S=S[[x]], alternative="greater", alpha=alpha), 
        extra_args
      ))
      return(discoveries(tmp))
    }
  }else{
    
    G <- flip::t2p(G, obs.only = FALSE)
    
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
  rownames(out) <- names(S)
  
  return(out)
}

