\name{sumSome-package}
\alias{sumSome-package}
\alias{sumSome}
\docType{package}
\title{
  True Discovery Guarantee by Sum-Based Tests
}
\description{
  It provides true discovery guarantees, using sum-based global statistics (sum of t-scores, p-value combinations, etc.). As main features, it produces permutation-based simultaneous lower confidence bounds for the proportion of active voxels in clusters for fMRI data, and for the proportion of differentially expressed genes in pathways for gene expression data.
}
\author{
Anna Vesely and Xu Chen.

Maintainer: Anna Vesely <anna.vesely@uni-bremen.de>
}
\references{
  Goeman J. J. and Solari A. (2011). Multiple testing for exploratory research. Statistical Science, doi: 10.1214/1-STS356.

  Tian J., Chen X., Katsevich E., Goeman J. J. and Ramdas A. (2022). Large-scale simultaneous inference under dependence. Scandinavian Journal of Statistics, doi: 10.1111/sjos.12614.

  Vesely A., Finos L., and Goeman J. J. (2023). Permutation-based true discovery guarantee by sum tests. Journal of the Royal Statistical Society, Series B (Statistical Methodology), doi: 10.1093/jrsssb/qkad019.
}
\keyword{ package }
\seealso{
  
  fMRI cluster analysis: \code{\link{brainScores}}, \code{\link{brainPvals}}, \code{\link{brainClusters}}, \code{\link{brainAnalysis}}
  
  Gene expression pathway analysis: \code{\link{geneScores}}, \code{\link{genePvals}}, \code{\link{geneAnalysis}}
  
  General setting: \code{\link{sumStats}} and \code{\link{sumPvals}} (permutations), \code{\link{sumStatsPar}} and \code{\link{sumPvalsPar}} (parametric)
}

