% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sumTest.R
\name{sumTest}
\alias{sumTest}
\title{True Discovery Guarantee}
\usage{
sumTest(G, S, alpha, truncFrom, truncTo, nMax)
}
\arguments{
\item{G}{numeric matrix of statistics, where columns correspond to variables, and rows to data transformations (e.g. permutations).
The first transformation is the identity. Extreme values are the greatest.}

\item{S}{vector of indices for the variables of interest.}

\item{alpha}{significance level.}

\item{truncFrom}{truncation parameter: values smaller than \code{truncFrom} are truncated.
If \code{NULL}, statistics are not truncated.}

\item{truncTo}{truncation parameter: truncated values are set to \code{truncTo}.
If \code{NULL}, statistics are not truncated.}

\item{nMax}{maximum number of iterations.}
}
\value{
\code{sumTest} returns an object of class \code{sumSome}, containing
\itemize{
\item \code{total}: total number of variables (columns in \code{G})
\item \code{size}: size of \code{S}
\item \code{alpha}: significance level
\item \code{TD}: lower (1-\code{alpha})-confidence bound for the number of true discoveries in \code{S}
\item \code{maxTD}: maximum value of \code{TD} that could be found under convergence of the algorithm
\item \code{iterations}: number of iterations of the algorithm
}
}
\description{
Internal function.
It determines confidence bounds for the number of true discoveries, the true discovery proportion
and the false discovery proportion within a set of interest.
The bounds are simultaneous over all sets, and remain valid under post-hoc selection.
}
\details{
Truncation parameters should be such that \code{truncTo} is not greater than \code{truncFrom}.

The significance level \code{alpha} should be in the interval [1/\code{B}, 1), where
\code{B} is the number of data transformations (rows in \code{G}).
}
\author{
Anna Vesely.
}
