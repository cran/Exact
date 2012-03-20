\name{Exact}
\alias{Exact}


\docType{package}
\title{Unconditional Exact Tests for 2x2 Tables}
\description{Uses Barnard's or Boschloo's test to compute the p-value of a given 2x2 table or calculate the power
for known proportions and sample sizes.}

\details{
\tabular{ll}{
Package: \tab Exact\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-16\cr
License: \tab GPL-2 \cr
}
Barnard's or Boschloo's test can be used to test the independence of rows and columns for a 2x2 table.
Unconditional tests are usually more powerful than conditional tests, but require more computation time.
P-values can be computed for 2x2 tables under the binomial and multinomial models with the option of 
using different statistics and using an interval method.  Power calculations can be computed for known sample
sizes and proportions. It has been shown that Barnard's test is uniformly more powerful than Fisher's test for
balanced designs.  There are several different methods for finding as or more extreme tables than the one
observed.  Suissa and Shuster suggested using a Z-pooled statistic.  Boschloo recommended using the p-value for 
Fisher's test as the test statistic.  This became known as Boschloo's test, which is always uniformly more powerful than
Fisher's test.  In addition, Berger and Boos proposed considering only values of the nuisance parameter that are in a 
constructed confidence interval.  They showed that for many cases, the interval approach often yield more powerful tests.
All of these tests can be done for both models, where the appropriate model is based on the design.
}
\author{
Peter Calhoun

Maintainer: Peter Calhoun <calhoun.peter@gmail.com>
}
\references{
Barnard, G.A. (1945) A new test for 2x2 tables. Nature, 156:177.

Barnard, G.A. (1947) Significance tests for 2x2 tables. Biometrika, 34:123-138.

Berger, R. and Boos D. (1994) P values maximized over a confidence set for the nuisance parameter. Journal of the American Statistical Association 89, 1012-1016.

Berger, R. (1994) Power comparison of exact unconditional tests for comparing two binomial proportions. Institute of Statistics Mimeo Series No. 2266.

Berger, R. (1996) More powerful tests from confidence interval p values. American Statistician, 50, 314-318.

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table
when Testing the Equality of Two Probabilities. Statistica Neerlandica, 24, 1-35.

Cardillo, G. (2009) MyBarnard: a very compact routine for Barnard's exact test on 2x2 matrix.
\url{http://www.mathworks.com/matlabcentral/fileexchange/25760}

Galili, T. and Calhoun, P. (2010)

\url{http://www.r-statistics.com/2010/02/barnards-exact-test-a-powerful-alternative-for-fishers-exact-test-implemented-in-r/}.

Suissa, S. and Shuster, J. J. (1985), Exact Unconditional Sample Sizes for the 2x2 Binomial Trial,
Journal of the Royal Statistical Society, Ser. A, 148, 317-327.

Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez, L. Rodriguez-Cardozo N.A. RamosDelgado and R. Garcia-Sanchez. (2004).
Barnardextest:Barnard's Exact Probability Test. A MATLAB file. [WWW document]
\url{http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=6198}.

}
