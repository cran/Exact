\name{exact.test}
\alias{exact.test}
\alias{binomialCode}
\alias{multinomialCode}
\title{Unconditional Exact Tests for 2x2 Tables with Independent Samples}
\description{Calculates Barnard's or Boschloo's unconditional exact test for binomial or multinomial models with independent samples}
\usage{
exact.test(data, alternative = c("two.sided", "less", "greater"), npNumbers = 100,
    np.interval = FALSE, beta = 0.001,
     method = c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm"),
    model = c("Binomial", "Multinomial"), tsmethod = c("square", "central"),
    conf.int = FALSE, conf.level = 0.95,
    cond.row = TRUE, to.plot = TRUE, ref.pvalue = TRUE,
    delta = 0, reject.alpha = NULL, useStoredCSM = TRUE)
}
\arguments{
  \item{data}{A two dimensional contingency table in matrix form}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding the more extreme tables:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", "Boschloo", or "CSM".
CSM tests cannot be calculated for multinomial models}
  \item{model}{The model being used: must be either "Binomial" or "Multinomial"}
  \item{tsmethod}{
Indicates two-sided method: must be either "square" or "central".  Only used if model="Binomial"}
  \item{conf.int}{
Logical: Indicates if a confidence interval on the difference in proportion should be computed.
Only used if model="Binomial"}
  \item{conf.level}{Number: Confidence level of interval on difference in proportion.  Only used if conf.int=TRUE}
  \item{cond.row}{Logical: Indicates if row margins are fixed in the binomial models.  Only used if model="Binomial"}
  \item{to.plot}{Logical: Indicates if plot of p-value vs. nuisance parameter should be generated.  Only used if model="Binomial"}
  \item{ref.pvalue}{Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected.  Only used if model="Binomial"}
  \item{delta}{Number: null hypothesis of the difference in proportion.  Only used if model="Binomial"}
  \item{reject.alpha}{Number: significance level for exact test.  Optional and primarily used to speed up \code{\link{exact.reject.region}} function}
  \item{useStoredCSM}{Logical: uses stored CSM ordering matrix.  Only used if method="csm"}
}

\details{
Unconditional exact tests (i.e., Barnard's test) can be performed for binomial or multinomial models with independent samples. The binomial model assumes the row or column margins (but not both) are known in advance, while the multinomial model assumes only the total sample size is known beforehand.  For the binomial model, the user needs to specify which margin is fixed (default is rows). Conditional tests (e.g., Fisher's exact test) have both row and column margins fixed, but this is a very uncommon design.  For paired samples, see \code{\link{paired.exact.test}}.

For the binomial model, the null hypothesis is the difference of proportion is equal to 0. Under the null hypothesis, the probability of a 2x2 table is the product of two binomials. The p-value is calculated by maximizing a nuisance parameter and summing the as or more extreme tables. The \code{method} parameter specifies the method to determine the more extreme tables (see references for more details):

\itemize{

  \item Z-pooled (or Score) - Uses the test statistic from a Z-test using a pooled proportion

  \item Z-unpooled - Uses the test statistic from a Z-test without using the pooled proportion
  
  \item Santner and Snell - Uses the difference in proportion
  
  \item Boschloo - Uses the p-value from Fisher's exact test
  
  \item CSM - Starts with the most extreme table and sequentially adds more extreme tables based on the smallest p-value (calculated by maximizing the probability of a 2x2 table).  This is Barnard's original method
}

There is some disagreement on which method to use.  Barnard's CSM (Convexity, Symmetry, and Minimization) test is often the most powerful test, but is much more computationally intensive.  This test is recommended by Mato and Andres and the author of this R package when computationally feasible.  Suissa and Shuster suggested using a Z-pooled statistic, which is uniformly more powerful than Fisher's test for balanced designs.  Boschloo recommended using the p-value for Fisher's test as the test statistic.  This method became known as Boschloo's test, and it is always uniformly more powerful than Fisher's test.  Many researchers argue that Fisher's exact test should \strong{never} be used to analyze 2x2 tables (except in the rare instance both margins are fixed).

Once the more extreme tables are determined, the p-value is calculated by maximizing over the common success probability -- a nuisance parameter.  The p-value computation has many local maxima and can be computationally intensive.  The code performs an exhaustive search by considering many values of the nuisance parameter from 0 to 1, represented by \code{npNumbers}.  If \code{ref.pvalue = TRUE}, then the code will also use the \code{optimise} function near the nuisance parameter to refine the p-value.  Increasing \code{npNumbers} and using \code{ref.pvalue} ensures the p-value is correctly calculated at the expense of slightly more computation time.

Another approach, proposed by Berger and Boos, is to calculate the Clopper-Pearson confidence interval of the nuisance parameter (represented by \code{np.interval}) and only maximize the p-value function for nuisance parameters within the confidence interval; this approach adds a small penalty to the p-value to control for the type 1 error rate (cannot be used with CSM).

The tests can also be implemented for non-inferiority hypotheses by changing the \code{delta} for binomial models. A confidence interval for the difference in proportion can be constructed by determining the smallest delta such that all greater deltas are significant (essentially the delta that crosses from non-significant to significant, but since the p-value is non-monotonic as a function of delta, this code takes the supremum).  For details regarding calculation, please see Chan (2003).  We note the "z-pooled" method uses a delta-projected Z-statistic (aka Score) and uses a constrained MLE of the success proportions, while "boschloo" method ignores the delta and uses the same ordering procedure.  

For two-sided tests, the code will either sum the probabilities for both sides of the table if \code{tsmethod = "square"} (Default, same approach as \code{fisher.test}) or construct a one-sided test and double the p-value if \code{tsmethod = "central"}.  Generally, the "square" procedure is more powerful and conventional, although there are some advantages with using the "central" procedure.  Mainly, Boschloo's test cannot order tables (i.e., use Fisher's p-value) for a two-sided alternative with non-zero delta.  Thus, to calculate a two-sided confidence interval with Boschloo's test, one has to resort to using the "central" approach.  For other tests, there is an equivalent statistic based on delta, and the two-sided p-value is determined by either the Agresti-Min interval (\code{tsmethod = "square"}) or Chan-Zhang interval (\code{tsmethod = "central"}).

The CSM test is computationally intensive due to iteratively maximizing the p-value calculation to order the tables.  The CSM ordering matrix has been stored for all possible sample sizes less than or equal to 100 (i.e., max(n1,n2)<=100).  Thus, using the \code{useStoredCSM = TRUE} can greatly improve computation time.  However, the stored ordering matrix was computed with \code{npNumbers=100} and it is possible that the ordering matrix was not optimal for larger \code{npNumbers}.  Increasing \code{npNumbers} and setting \code{useStoredCSM = FALSE} ensures the p-value is correctly calculated at the expense of significantly greater computation time.  The stored ordering matrix is not used in the calculation of confidence intervals or non-inferiority tests, so CSM can still be very computationally intensive.

The above description applies to the binomial model.  The multinomial model is similar except there are two nuisance parameters.  The CSM test has not been developed for multinomial models.  Improvements to the code have focused on the binomial model, so multinomial models takes substantially longer.
}

\value{
A list with class "htest" containing the following components:

\item{p.value}{The p-value of the test}
\item{statistic}{The observed test statistic to determine more extreme tables}
\item{estimate}{An estimate of the difference in proportion}
\item{null.value}{The difference in proportion under the null}
\item{conf.int}{A confidence interval for the difference in proportion.  Only present if \code{conf.int = TRUE}}
\item{alternative}{A character string describing the alternative hypothesis}
\item{model}{A character string describing the model design ("Binomial" or "Multinomial")}
\item{method}{A character string describing the method to determine more extreme tables}
\item{tsmethod}{A character string describing the method to implement two-sided tests.  Only present if \code{model = "binomial"}}
\item{np}{The nuisance parameter that maximizes the p-value.  For multinomial models, both nuisance parameters are given}
\item{np.range}{The range of nuisance parameters considered.  For multinomial models, both nuisance parameter ranges are given}
\item{data.name}{A character string giving the names of the data}
\item{parameter}{The sample sizes}
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
\url{https://www.mathworks.com/matlabcentral/fileexchange/25760-mybarnard}

Mato, S. and Andres, M. (1997), Simplifying the calculation of the P-value for Barnard's test and its derivatives.  \emph{Statistics and Computing}, \bold{7}, 137--143

Mehrotra, D., Chan, I., Berger, R. (2003), A Cautionary Note on Exact Unconditional Inference for a Difference Between Two Independent Binomial Proportions.  \emph{Biometrics}, \bold{59}, 441--450

Ruxton, G. D. and Neuhauser M (2010), Good practice in testing for an association in contingency tables. \emph{Behavioral Ecology and Sociobiology}, \bold{64}, 1505--1513

Suissa, S. and Shuster, J. J. (1985), Exact Unconditional Sample Sizes for the 2x2 Binomial Trial, \emph{Journal of the Royal Statistical Society, Ser. A}, \bold{148}, 317--327
}
\author{Peter Calhoun}
\note{
CSM test and multinomial models are much more computationally intensive.  I have also spent a greater amount of time improving the computation speed and adding functionality to the binomial model. Increasing the number of nuisance parameters considered and refining the p-value will increase the computation time, but more likely to ensure accurate calculations. Performing confidence intervals also greatly increases computation time.

This code was influenced by the FORTRAN program located at \url{https://www4.stat.ncsu.edu/~boos/exact/}
}
\section{Warning}{Multinomial models and CSM tests with confidence intervals may take a very long time, even for small sample sizes.}

\seealso{
\code{fisher.test} and \pkg{exact2x2}
}
\examples{
data <- matrix(c(7, 8, 12, 3), 2, 2, byrow=TRUE)
exact.test(data, alternative="two.sided", method="Z-pooled",
           conf.int=TRUE, conf.level=0.95)
exact.test(data, alternative="two.sided", method="Boschloo",
           np.interval=TRUE, beta=0.001, tsmethod="central")

\dontrun{
# Ensure that the ExactData R package is available before running the CSM test.
if (requireNamespace("ExactData", quietly = TRUE)) {
# Example from Barnard's (1947) appendix:
data <- matrix(c(4, 0, 3, 7), 2, 2,
               dimnames=list(c("Box 1","Box 2"), c("Defective","Not Defective")))
exact.test(data, method="CSM", alternative="two.sided")
}
}
}
\keyword{htest}
\concept{Barnard}
\concept{Boschloo}
\concept{Unconditional}
\concept{Exact Tests}