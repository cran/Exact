\name{Exact}
\alias{Exact}
\alias{binom.CI}
\alias{boschloo_TX}
\alias{csm_Tbls}
\alias{csmApprox_TX}
\alias{csmMod_Tbls}
\alias{fisher.2x2}
\alias{maxPvalue}
\alias{santner_TX}
\alias{zpooled_TX}
\alias{zunpooled_TX}

\docType{package}
\title{Unconditional Exact Tests for 2x2 Tables}
\description{
This package implements the \code{\link{exact.test}} function to perform unconditional exact tests.  This package also includes the \code{\link{power.exact.test}} function to calculate the power to detect a significant difference using unconditional exact tests.
}


\details{

\tabular{ll}{
Package: \tab Exact\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2016-10-21\cr
License: \tab GPL-2 \cr
}

Unconditional exact tests can be performed to test the independence of rows and columns in a 2x2 table.  Unconditional tests (such as Barnard's and Boschloo's exact tests) are more powerful alternatives than conditional tests (such as Fisher's exact test).  P-values can be computed for 2x2 tables under the binomial and multinomial models using various test statistics to find the 'as or more extreme' tables.  The multinomial model assumes only the total sample size is known in advance (common in cross-sectional studies), the binomial model assumes only the row or column margins (but not both) are fixed (most common situation; occurs in case-control studies), and the hypergeometric model assumes both row and column margins are fixed (very unlikely situation; only design where Fisher's test should be performed).

The details of the test statistics are given in the \code{\link{exact.test}} function description.  Suissa and Shuster suggested using a Z-pooled statistic, which they found to be uniformly more powerful than Fisher's test for balanced designs.  Boschloo recommended using the p-value for Fisher's test as the test statistic.  This method became known as Boschloo's test, and it is always uniformly more powerful than Fisher's test.  Mato and Andres suggested using Barnard's CSM test. Additionally, Berger and Boos proposed considering only values of the nuisance parameter that are in a constructed confidence interval.  The interval approach often yields more powerful tests.  While there is still a disagreement on which test statistic to use, most researchers agree that Fisher's exact test should \strong{not} be used to analyze 2x2 tables.  All of these tests can be computed in this package.
}

\note{
Thoughout the years I have received help while creating this package. 
Special thanks goes to Philo Calhoun, Tal Galili, Kamil Erguler, Roger Berger, Karl Hufthammer, and the R community.
}

\author{
\packageAuthor{Exact}

Maintainer: \packageMaintainer{Exact}
}
\references{
Barnard, G.A. (1945) A new test for 2x2 tables. \emph{Nature}, \bold{156}, 177

Barnard, G.A. (1947) Significance tests for 2x2 tables. \emph{Biometrika}, \bold{34}, 123-138

Berger, R. and Boos D. (1994) P values maximized over a confidence set for the nuisance parameter. \emph{Journal of the American Statistical Association}, \bold{89}, 1012-1016

Berger, R. (1994) Power comparison of exact unconditional tests for comparing two binomial proportions. \emph{Institute of Statistics Mimeo Series No. 2266}

Berger, R. (1996) More powerful tests from confidence interval p values. \emph{American Statistician}, \bold{50}, 314-318

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table when Testing the Equality of Two Probabilities. \emph{Statistica Neerlandica}, \bold{24}, 1-35

Cardillo, G. (2009) MyBarnard: a very compact routine for Barnard's exact test on 2x2 matrix.
\url{http://www.mathworks.com/matlabcentral/fileexchange/25760-mybarnard}

Mato, S. and Andres, M. (1997), Simplifying the calculation of the P-value for Barnard's test and its derivatives.  \emph{Statistics and Computing}, \bold{7}, 137-143

Mehrotra, D., Chan, I., Berger, R. (2003), A Cautionary Note on Exact Unconditional Inference for a Difference Between Two Independent Binomial Proportions.  \emph{Biometrics}, \bold{59}, 441-450

Ruxton, G. D. and Neuhauser M (2010), Good practice in testing for an association in contingency tables. \emph{Behavioral Ecology and Sociobiology}, \bold{64}, 1505-1513

Suissa, S. and Shuster, J. J. (1985), Exact Unconditional Sample Sizes for the 2x2 Binomial Trial, \emph{Journal of the Royal Statistical Society, Ser. A}, \bold{148}, 317-327
}

\keyword{Exact Tests}
