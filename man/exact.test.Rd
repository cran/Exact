\name{exact.test}
\alias{exact.test}
\alias{binomialCode}
\alias{multinomialCode}
\title{Unconditional Exact Tests for 2x2 tables}
\description{Calculates Barnard's or Boschloo's unconditional exact test for binomial or multinomial models}
\usage{
exact.test(data, alternative = c("two.sided", "less", "greater"), npNumbers = 100,
    np.interval = FALSE, beta = 0.001,
    method = c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
               "csm", "csm approximate"),
    model = c("Binomial", "Multinomial"), conf.int = FALSE, conf.level = 0.95,
    precision = 0.001, cond.row = TRUE, to.plot = TRUE, ref.pvalue = TRUE,
    delta = 0, reject.alpha = NULL)
}

\arguments{
  \item{data}{A two dimensional contingency table in matrix form}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding the more extreme tables:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", "Boschloo", "CSM", or "CSM approximate".
CSM tests cannot be calculated for multinomial models
}
  \item{model}{The model being used: must be either "Binomial" or "Multinomial"}
  \item{conf.int}{
Logical: Indicates if a confidence interval on the difference in proportion should be computed.
Only used if model="Binomial" and method is "z-pooled" or "csm"
}
  \item{conf.level}{Number: Confidence level of interval on difference in proportion.  Only used if conf.int=TRUE}
  \item{precision}{
Number: Precision of confidence interval on difference in proportion.  Smaller number represents more digits in confidence interval.  Only used if conf.int=TRUE
}
  \item{cond.row}{Logical: Indicates if row margins are fixed in the binomial models.  Only used if model="Binomial"}
  \item{to.plot}{Logical: Indicates if plot of p-value vs. nuisance parameter should be generated.  Only used if model="Binomial"}
  \item{ref.pvalue}{Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected.  Only used if model="Binomial"}
  \item{delta}{Number: null hypothesis of the difference in proportion.  Only used if model="Binomial" and method is "z-pooled" or "csm"}
  \item{reject.alpha}{Number: significance level for exact test.  Optional and primarily used to speed up \code{\link{exact.reject.region}} function}
}
\details{
Unconditional exact tests can be performed for binomial or multinomial models. The binomial model assumes the row or column margins (but not both) are known in advance, while the multinomial model assumes only the total sample size is known beforehand.  For the binomial model, the user needs to specify which margin is fixed (default is rows). Conditional tests (e.g., Fisher's exact test) have both row and column margins fixed, but this is a very uncommon design.

For the binomial model, the null hypothesis is the difference of proportion is equal to 0. Under the null hypothesis, the probability of a 2x2 table is the product of two binomials. The p-value is calculated by maximizing a nuisance parameter and summing the as or more extreme tables. The \code{method} parameter specifies the method to determine the more extreme tables (see references for more details):

\itemize{

  \item Z-pooled (or Score) - Uses the test statistic from a Z-test using a pooled proportion

  \item Z-unpooled (or Wald) - Uses the test statistic from a Z-test without using the pooled proportion
  
  \item Santner and Snell - Uses the difference in proportion
  
  \item Boschloo - Uses the p-value from Fisher's exact test
  
  \item CSM - Starts with the most extreme table and sequentially adds more extreme tables based on the smallest p-value (calculated by maximizing the probability of a 2x2 table).  This is Barnard's original method
  
  \item CSM approximate - Similar to the CSM test, but does not update the maximum p-value calculation when each table is added
}

There is some disagreement on which method to use.  Suissa and Shuster suggested using a Z-pooled statistic, which is uniformly more powerful than Fisher's test for balanced designs.  Boschloo recommended using the p-value for Fisher's test as the test statistic.  This method became known as Boschloo's test, and it is always uniformly more powerful than Fisher's test.  Mato and Andres suggested using Barnard's CSM test.  This is the recommended method by the author of this R package, but is much more computationally intensive.  Regardless of the method chosen, most researchers agree that Fisher's exact test should \strong{not} be used to analyze 2x2 tables.

Once the more extreme tables are determined, the p-value is calculated the maximizing over the common success probability -- a nuisance parameter.  The p-value computation has many local maxima and can be computationally intensive.  The code performs an exhaustive search by considering many values of the nuisance parameter from 0 to 1, represented by \code{npNumbers}.  Another approach, proposed by Berger and Boos, is to calculate the Clopper-Pearson confidence interval of the nuisance parameter (represented by \code{np.interval}) and only maximizing the p-value function for nuisance parameters within the confidence interval; this approach adds a small penalty to the p-value to control for the type 1 error rate (cannot be used with CSM).  If \code{ref.pvalue = TRUE}, then the code will also use the \code{optimise} function near the nuisance parameter to refine the p-value.  Increasing \code{npNumbers} and using \code{ref.pvalue} ensures the p-value is correctly calculated at the expense of slightly more computation time.

The Z-pooled and CSM test has been extended to test the difference in proportion is equal to \code{delta} (not only 0).  Thus, the Z-pooled and CSM test can be used to test noninferiority designs, and be used to construct confidence intervals for the difference in proportion.

There are many ways to define the two-sided p-value; this code uses the \code{fisher.test} approach by summing the probabilities for both sides of the table.  Confidence intervals of the difference in proportion (Z-pooled and CSM test only) use the Agresti-Min interval approach for two-sided tests and the Chan-Zhang interval approach for one-sided tests.  The \code{precision} represents the number of decimal points for the confidence interval.

The above description applies to the binomial model.  The multinomial model is similar except there are two nuisance parameters.  The CSM test has not been developed for multinomial models nor have noninferiority tests been developed.  Improvements to the code have focused on the binomial model, so multinomial models takes substantially longer.
}

\value{
A list with class "htest" containing the following components:

\item{p.value}{The computed p-value}
\item{test.statistic}{The observed test statistic}
\item{estimate}{An estimate of the parameter tested}
\item{alternative}{A character string describing the alternative hypothesis}
\item{model}{A character string describing the model design ("Binomial" or "Multinomial")}
\item{method}{A character string describing the method to determine more extreme tables}
\item{np}{The nuisance parameter that maximizes the p-value.  For multinomial models, both nuisance parameters are given}
\item{np.range}{The range of nuisance parameters considered.  For multinomial models, both nuisance parameter ranges are given}
\item{data.name}{A character string giving the names of the data}
}

\references{
Agresti, A. and Min, Y. (2001) On small-sample confidence intervals for parameters in discrete distributions. \emph{Biometrics}, \bold{57}, 963--971

Barnard, G.A. (1945) A new test for 2x2 tables. \emph{Nature}, \bold{156}, 177

Barnard, G.A. (1947) Significance tests for 2x2 tables. \emph{Biometrika}, \bold{34}, 123--138

Berger, R. and Boos D. (1994) P values maximized over a confidence set for the nuisance parameter. \emph{Journal of the American Statistical Association}, \bold{89}, 1012--1016

Berger, R. (1994) Power comparison of exact unconditional tests for comparing two binomial proportions. \emph{Institute of Statistics Mimeo Series No. 2266}

Berger, R. (1996) More powerful tests from confidence interval p values. \emph{American Statistician}, \bold{50}, 314--318

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table when Testing the Equality of Two Probabilities. \emph{Statistica Neerlandica}, \bold{24}, 1--35

Chan, I. (2003), Proving non-inferiority or equivalence of two treatments with dichotomous endpoints using exact methods, \emph{Statistical Methods in Medical Research}, \bold{12}, 37--58

Cardillo, G. (2009) MyBarnard: a very compact routine for Barnard's exact test on 2x2 matrix.
\url{http://www.mathworks.com/matlabcentral/fileexchange/25760-mybarnard}

Mato, S. and Andres, M. (1997), Simplifying the calculation of the P-value for Barnard's test and its derivatives.  \emph{Statistics and Computing}, \bold{7}, 137--143

Mehrotra, D., Chan, I., Berger, R. (2003), A Cautionary Note on Exact Unconditional Inference for a Difference Between Two Independent Binomial Proportions.  \emph{Biometrics}, \bold{59}, 441--450

Ruxton, G. D. and Neuhauser M (2010), Good practice in testing for an association in contingency tables. \emph{Behavioral Ecology and Sociobiology}, \bold{64}, 1505--1513

Suissa, S. and Shuster, J. J. (1985), Exact Unconditional Sample Sizes for the 2x2 Binomial Trial, \emph{Journal of the Royal Statistical Society, Ser. A}, \bold{148}, 317--327
}
\author{Peter Calhoun}
\note{
CSM test and multinomial models are much more computationally intensive.  I have also spent a greater amount of time making the computations for the binomial models more efficient; future work will be devoted to improving the multinomial models. Increasing the number of nuisance parameters considered and refining the p-value will increase the computation time. Performing confidence intervals also greatly increases computation time.

This code was influenced by the FORTRAN program located at \url{http://www4.stat.ncsu.edu/~boos/exact/}
}
\section{Warning}{Multinomial models and CSM tests may take a very long time, even for small sample sizes.}

\seealso{
\code{fisher.test} and \pkg{exact2x2}
}
\examples{
data <- matrix(c(7, 8, 12, 3), 2, 2, byrow=TRUE)
exact.test(data, alternative="less",to.plot=TRUE)
exact.test(data, method="Z-pooled", alternative="two.sided", conf.int=TRUE, conf.level=0.95,
           precision=0.0001, npNumbers=100)
exact.test(data, method="Boschloo", alternative="two.sided", np.interval=TRUE, beta=0.001,
           npNumbers=100, to.plot=FALSE)

# Example from Barnard's (1947) appendix:
data <- matrix(c(4, 0, 3, 7), 2, 2,
               dimnames=list(c("Box 1","Box 2"), c("Defective","Not Defective")))
exact.test(data, method="CSM", alternative="two.sided", conf.int=TRUE)

data <- matrix(c(6, 8, 4, 3), 2, 2, byrow=TRUE)
exact.test(data, model="Multinomial", alternative="less", method="Z-pooled")
}
\keyword{Barnard}
\keyword{Boschloo}
\keyword{Unconditional}
\keyword{Exact Tests}

