\name{sumSome-package}
\alias{sumSome-package}
\alias{sumSome}
\docType{package}
\title{
  True Discovery Guarantee by Sum-Based Tests
}
\description{
  It provides permutation-based true discovery guarantees, using sum-based global statistics 
  (sum of t-scores, p-value combinations, etc.). As a main feature, it produces simultaneous lower confidence bounds
  for the number of active voxels in different clusters for fMRI cluster analysis.
}
\author{
Anna Vesely.

Maintainer: Anna Vesely <anna.vesely@phd.unipd.it>
}
\references{

  Goeman, J. J., and Solari, A. (2011). Multiple testing for exploratory research. Statistical Science 26 (4) 584-597.
  
  Vesely, A., Finos, L., and Goeman, J. J. (2020). Permutation-based true discovery guarantee by sum tests. Pre-print arXiv:2102.11759.
}
\keyword{ package }
\seealso{
  True discovery guarantee: \code{\link{sumStats}}, \code{\link{sumPvals}}
  
  Permutation statistics for brain imaging: \code{\link{brainScores}}, \code{\link{brainPvals}}
  
  True discovery guarantee for cluster analysis: \code{\link{clusterAnalysis}}
}

