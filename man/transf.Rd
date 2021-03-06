% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transf.R
\name{transf}
\alias{transf}
\title{Transformation of Statistics}
\usage{
transf(G, truncFrom, truncTo, option, r)
}
\arguments{
\item{G}{numeric matrix of statistics.}

\item{truncFrom}{truncation parameter: values less extreme than \code{truncFrom} are truncated.
If \code{NULL}, statistics are not truncated.}

\item{truncTo}{truncation parameter: truncated values are set to \code{truncTo}.
If \code{NULL}, statistics are not truncated.}

\item{option}{direction of the alternative hypothesis (\code{greater}, \code{lower}, \code{two.sided}),
or transformation (\code{squares} for generic statistics,
and \code{edgington}, \code{fisher}, \code{pearson}, \code{liptak}, \code{cauchy}, \code{vovk.wang} for p-values).}

\item{r}{parameter for Vovk and Wang's p-value combination.}
}
\value{
\code{transf} returns a list containing the truncated and transformed matrix \code{G},
and the transformed truncation parameters \code{truncFrom} and \code{truncTo}.
}
\description{
Internal function.
It truncates and transforms a matrix of statistics.
}
\details{
Transformations are defined so that the most extreme values of the new statistics are always the greatest.
A generic statistic \code{x} is transformed as following.
\itemize{
\item greater: \code{x}
\item lower: \code{-x}
\item two-sided: \code{|x|}
\item squares: \code{x^2}
\item Edgington: \code{-x}
\item Fisher: \code{-log(x)}
\item Pearson: \code{log(1-x)}
\item Liptak: \code{-qnorm(x)}
\item Cauchy: \code{tan(0.5 - x)/x}
\item Vovk and Wang: \code{- sign(r)x^r}
}
An error message is returned if the transformation produces infinite values.

Truncation parameters should be such that \code{truncTo} is not more extreme than \code{truncFrom}.
As Pearson's and Liptak's transformations produce infinite values in 1, for such methods
\code{truncTo} should be strictly smaller than 1.
}
\author{
Anna Vesely.
}
\keyword{Internal}
