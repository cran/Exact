\name{Exact-package}
\alias{Exact}
\alias{binom.CI}
\alias{checkPairedParam}
\alias{checkParam}
\alias{chisq_TX}
\alias{confInt}
\alias{confIntPaired}
\alias{convertMethod}
\alias{csmPairedTemp}
\alias{csmPairedTemp2sidedDelta}
\alias{csmTemp}
\alias{csmTemp2sidedDelta}
\alias{dbinomCalc}
\alias{discordant.CI}
\alias{fisher.2x2}
\alias{lowerMatVal}
\alias{maxPvalue}
\alias{maxPvalueLookup}
\alias{maxPvaluePaired}
\alias{maxPvaluePairedLookup}
\alias{maxPvalueMultinom}
\alias{mcnemar.2x2}
\alias{mcnemar_TX}
\alias{methodText}
\alias{moreExtreme}
\alias{moreExtremeCSM}
\alias{moreExtremePaired}
\alias{moreExtremePairedCSM}
\alias{multinom}
\alias{pairedMethodText}
\alias{santner_TX}
\alias{searchExtreme}
\alias{searchExtremePaired}
\alias{trinom}
\alias{trinomCalc}
\alias{updateMat}
\alias{zpooled_TX}
\alias{zunpooled_TX}

\docType{package}
\title{Unconditional Exact Tests R Package}
\description{
This package performs unconditional exact tests using \code{\link{exact.test}} for independent samples or \code{\link{paired.exact.test}} function for paired samples.  The unconditional exact tests for independent samples are referred to as Barnard's (1945, 1947) test and also extended to test two paired proportions.  This package also includes the \code{\link{power.exact.test}} and \code{\link{power.paired.test}} functions to calculate the power of various tests and the \code{\link{exact.reject.region}} and \code{\link{paired.reject.region}} functions to determine the rejection region of a test.
}
\details{
Unconditional exact tests are a more powerful alternative than conditional exact tests.  This package can compute p-values, confidence intervals, and power calculations for various tests.  Details of the tests are given in the \code{\link{exact.test}} documentation for independent samples and \code{\link{paired.exact.test}} documentation for paired samples.
}
\author{
\packageAuthor{Exact}

Maintainer: \packageMaintainer{Exact}
}
\references{
Barnard, G.A. (1945) A new test for 2x2 tables. \emph{Nature}, \bold{156}, 177

Barnard, G.A. (1947) Significance tests for 2x2 tables. \emph{Biometrika}, \bold{34}, 123--138

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table when Testing the Equality of Two Probabilities. \emph{Statistica Neerlandica}, \bold{24}, 1--35

Berger, R.L. and Sidik, K. (2003) Exact unconditional tests for 2 x 2 matched-pairs design. \emph{Statistical Methods in Medical Research}, \bold{12}, 91--108

Hsueh, H., Liu, J., and Chen, J.J. (2001) Unconditional exact tests for equivalence or noninferiority for paired binary endpoints. \emph{Biometrics}, \bold{57}, 478--483
}
\note{
Thoughout the years I have received help while creating this package. 
Special thanks goes to Philo Calhoun, Tal Galili, Kamil Erguler, Roger Berger, Karl Hufthammer, and the R community.
}
\keyword{package}
\keyword{htest}
\concept{Barnard}
\concept{Boschloo}
\concept{McNemar}
\concept{Unconditional}
\concept{Exact Tests}