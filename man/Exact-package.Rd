\name{Exact-package}
\alias{Exact}
\alias{binom.CI}
\alias{chisq_TX}
\alias{confInt}
\alias{csmApprox_TX}
\alias{csmTemp}
\alias{csmTemp2sidedDelta}
\alias{dbinomCalc}
\alias{fisher.2x2}
\alias{maxPvalue}
\alias{maxPvalueLookup}
\alias{methodText}
\alias{moreExtreme}
\alias{moreExtremeCSM}
\alias{santner_TX}
\alias{searchExtreme}
\alias{zpooled_TX}
\alias{zunpooled_TX}

\docType{package}
\title{Unconditional Exact Tests R Package}
\description{
This package performs Barnard's (1945, 1947) or Boschloo's (1970) unconditional exact tests using the \code{\link{exact.test}} function.  This package also includes the \code{\link{power.exact.test}} function to calculate the power of various tests and the \code{\link{exact.reject.region}} to determine the rejection region of a test for a given sample size.
}
\details{
Unconditional exact tests are a more powerful alternative than conditional exact tests.  This package can compute p-values, confidence intervals, and power calculations for various tests.  Details of the tests are given in the \code{\link{exact.test}} documentation.
}
\references{
Barnard, G.A. (1945) A new test for 2x2 tables. \emph{Nature}, \bold{156}, 177

Barnard, G.A. (1947) Significance tests for 2x2 tables. \emph{Biometrika}, \bold{34}, 123--138

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table when Testing the Equality of Two Probabilities. \emph{Statistica Neerlandica}, \bold{24}, 1--35
}
\author{
\packageAuthor{Exact}

Maintainer: \packageMaintainer{Exact}
}
\note{
Thoughout the years I have received help while creating this package. 
Special thanks goes to Philo Calhoun, Tal Galili, Kamil Erguler, Roger Berger, Karl Hufthammer, and the R community.
}
\keyword{package}
\keyword{htest}
\concept{Exact Tests}