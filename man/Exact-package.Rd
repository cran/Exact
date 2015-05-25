\name{Exact}
\alias{Exact}
\alias{MultinomialCode}
\alias{BinomialCode}
\alias{fisher.2x2}

\docType{package}
\title{
Unconditional Exact Tests for 2x2 Tables
}
\description{This package implements the \code{exact.test} function to perform unconditional exact tests.  This package also
includes the \code{power.exact.test} function to calculate the power to detect a significant difference using unconditional exact
tests.}

\details{
\tabular{ll}{
Package: \tab Exact\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2015-05-25\cr
License: \tab GPL-2 \cr
}
Unconditional exact tests can be performed to test the independence of rows and columns in a 2x2 table.
Unconditional tests (such as Barnard's and Boschloo's exact tests) are more powerful alternatives than conditional tests (such as Fisher's exact test).
P-values can be computed for 2x2 tables under the binomial and multinomial models using various test statistics 
to find the 'as or more extreme' tables.  The multinomial model assumes only the total sample size is known in advance (common in cross-sectional studies),
the binomial model assumes only the row or column margins (but not both) are fixed (most common situation; occurs in case-control studies), and 
the hypergeometric model assumes both row and column margins are fixed (very unlikely situation; only design where Fisher's test should be performed).


The details of the test statistics are given in the \code{exact.test} function description.
Suissa and Shuster suggested using a Z-pooled statistic, which they found to be uniformly more powerful than Fisher's test 
for balanced designs.  Boschloo recommended using the p-value for Fisher's test as the test statistic.  This method became 
known as Boschloo's test, and it is always uniformly more powerful than Fisher's test.  Mato and Andres suggested using Barnard's CSM test.
Additionally, Berger and Boos proposed considering only values of the nuisance parameter that are in a constructed confidence interval.
The interval approach often yields more powerful tests.  While there is still a disagreement on which test statistic to use, most researchers agree
that Fisher's exact test should \emph{not} be used to analyze 2x2 tables.  All of these tests can be computed in this package.
}

\note{
Thoughout the years I have received help while creating this package. 
Special thanks goes to Philo Calhoun, Tal Galili, Kamil Erguler, Roger Berger, Karl Hufthammer, and the R community.
}

\author{
Peter Calhoun

Maintainer: Peter Calhoun <calhoun.peter@gmail.com>
}

\references{
Barnard, G.A. (1945) A new test for 2x2 tables. Nature, 156:177

Barnard, G.A. (1947) Significance tests for 2x2 tables. Biometrika, 34:123-138

Berger, R. and Boos D. (1994) P values maximized over a confidence set for the nuisance parameter. Journal of the American Statistical Association 89, 1012-1016

Berger, R. (1994) Power comparison of exact unconditional tests for comparing two binomial proportions. Institute of Statistics Mimeo Series No. 2266

Berger, R. (1996) More powerful tests from confidence interval p values. American Statistician, 50, 314-318

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table
when Testing the Equality of Two Probabilities. Statistica Neerlandica, 24, 1-35

Cardillo, G. (2009) MyBarnard: a very compact routine for Barnard's exact test on 2x2 matrix.
\url{http://www.mathworks.com/matlabcentral/fileexchange/25760-mybarnard}

Mato, S. and Andres, M. (1997), Simplifying the calculation of the P-value for Barnard's test and its derivatives.  Statistics and Computing, 7, 137-143

Mehrotra, D., Chan, I., Berger, R. (2003), A Cautionary Note on Exact Unconditional Inference for a Difference Between Two Independent Binomial Proportions.  Biometrics, 59, 441-450

Ruxton, G. D. and Neuhauser M (2010), Good practice in testing for an association in contingency tables. Behavioral Ecology and Sociobiology, 64, 1505-1513

Suissa, S. and Shuster, J. J. (1985), Exact Unconditional Sample Sizes for the 2x2 Binomial Trial,
Journal of the Royal Statistical Society, Ser. A, 148, 317-327

Trujillo-Ortiz, A., Hernandez-Walls, R., Castro-Perez, A., Rodriguez-Cardozo, L, RamosDelgado, N.A., Garcia-Sanchez, R. (2004).
Barnardextest:Barnard's Exact Probability Test. A MATLAB file. [WWW document]
\url{http://www.mathworks.com/matlabcentral/fileexchange/6198-barnardextest}

}
