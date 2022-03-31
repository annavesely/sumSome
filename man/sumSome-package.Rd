\name{sumSome-package}
\alias{sumSome-package}
\alias{sumSome}
\docType{package}
\title{
  True Discovery Guarantee by Sum-Based Tests
}
\description{
  It provides true discovery guarantees, using sum-based global statistics (sum of t-scores, p-value combinations, etc.). As main features, it produces simultaneous lower confidence bounds for the proportion of active voxels in clusters for fMRI data, and for the proportion of differentially expressed genes in pathways for gene expression data.
}
\author{
Anna Vesely and Xu Chen.

Maintainer: Anna Vesely <anna.vesely@unipd.it>
}
\references{

  Goeman, J. J., and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science 26 (4) 584-597.
  
  Hemerik, J. and Goeman, J. J. (2018). False discovery proportion estimation by permutations: confidence for significance analysis of microarrays. JRSS B, 80(1):137-155.
  
  Tian, J., Chen, X., Katsevich, E., Goeman, J. J., and Ramdas, A. (2021). Large-scale simultaneous inference under dependence. Pre-print arXiv:2102.11253.
  
  Vesely, A., Finos, L., and Goeman, J. J. (2021). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
}
\keyword{ package }
\seealso{
  
  fMRI cluster analysis: \code{\link{brainScores}}, \code{\link{brainPvals}}, \code{\link{brainClusters}}, \code{\link{brainAnalysis}}
  
  Gene expression pathway analysis: \code{\link{geneScores}}, \code{\link{genePvals}}, \code{\link{geneAnalysis}}
  
  General setting : \code{\link{sumStats}} and \code{\link{sumPvals}} (permutations), \code{\link{psumStats}} and \code{\link{psumPvals}} (parametric)
}

