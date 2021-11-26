\name{paired.exact.test}
\alias{paired.exact.test}
\alias{pairedCode}
\title{Unconditional Exact Tests for 2x2 Tables with Paired Samples}

\description{Calculates unconditional exact test for paired samples}
\usage{
paired.exact.test(data, alternative = c("two.sided", "less", "greater"), npNumbers = 100,
    np.interval = FALSE, beta = 0.001,
    method = c("uam", "ucm", "uamcc", "csm"),
    tsmethod = c("square", "central"),
    conf.int = FALSE, conf.level = 0.95, to.plot = TRUE,
    ref.pvalue = TRUE, delta = 0,
    reject.alpha = NULL, useStoredCSM = TRUE)
}
\arguments{
  \item{data}{A two dimensional contingency table in matrix form}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding the more extreme tables:
must be either "UAM", "UCM", "UAMCC", or "CSM"}
  \item{tsmethod}{A character string describing the method to implement two-sided tests.}
  \item{conf.int}{
Logical: Indicates if a confidence interval on the difference in proportion should be computed}
  \item{conf.level}{Number: Confidence level of interval on difference in proportion.  Only used if conf.int=TRUE}
  \item{to.plot}{Logical: Indicates if plot of p-value vs. nuisance parameter should be generated}
  \item{ref.pvalue}{Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected}
  \item{delta}{Number: null hypothesis of the difference in proportion}
  \item{reject.alpha}{Number: significance level for exact test.  Optional and primarily used to speed up \code{\link{paired.reject.region}} function}
  \item{useStoredCSM}{Logical: uses stored CSM ordering matrix.  Only used if method="csm"}
}
\details{
This function performs unconditional exact tests to compare two paired proportions.  The null hypothesis is the difference of two paired proportions is equal to 0. Under the null hypothesis, the probability of a 2x2 table is a trinomial distribution. The p-value is calculated by maximizing a nuisance parameter and summing the as or more extreme tables. The \code{method} parameter specifies the method to determine the more extreme tables (see references for more details):

\itemize{

  \item UAM (Unconditional Asymptotic McNemar) - Uses McNemar's Z-statistic

  \item UCM (Unconditional Conditional McNemar) - Uses McNemar's conditional p-value

  \item UAMCC (Unconditional Asyptotic McNemar with Continuity Correction) - Uses McNemar's Z-statistic with Continuity Correction
  
  \item CSM - Starts with the most extreme table and sequentially adds more extreme tables based on the smallest p-value (calculated by maximizing the probability of a 2x2 table).  This is extending Barnard's original method to test two paired proportions
}

There is little research comparing two paired proportions.  The author of this R package recommends using the CSM (Convexity, Symmetry, and Minimization) test as this test is often the most powerful, but is much more computationally intensive.

Once the more extreme tables are determined, the p-value is calculated by maximizing over the common discordant probability -- a nuisance parameter.  The p-value computation has many local maxima and can be computationally intensive.  The code performs an exhaustive search by considering many values of the nuisance parameter from 0 to 0.5, represented by \code{npNumbers}.  If \code{ref.pvalue = TRUE}, then the code will also use the \code{optimise} function near the nuisance parameter to refine the p-value.  Increasing \code{npNumbers} and using \code{ref.pvalue} ensures the p-value is correctly calculated at the expense of slightly more computation time.

Another approach, proposed by Berger and Sidik, is to calculate the Clopper-Pearson confidence interval of the nuisance parameter (represented by \code{np.interval}) and only maximize the p-value function for nuisance parameters within the confidence interval; this approach adds a small penalty to the p-value to control for the type 1 error rate (cannot be used with CSM).  If the CSM test is too computationally intensive, the author of this R package generally recommends using the UAM test with confidence interval approach.

The tests can also be implemented for non-inferiority hypotheses by changing the \code{delta}. A confidence interval for the difference in proportion can be constructed by determining the smallest delta such that all greater deltas are significant (essentially the delta that crosses from non-significant to significant, but since the p-value is non-monotonic as a function of delta, this code takes the supremum).  We note the "UAM" method uses a delta-projected Z-statistic, while "UCM" method ignores the delta and uses the same ordering procedure.

For two-sided tests, the code will either sum the probabilities for both sides of the table if \code{tsmethod = "square"} (Default) or construct a one-sided test and double the p-value if \code{tsmethod = "central"}.  The two methods give the same results when delta is zero.  For a non-zero delta, the "square" procedure is generally more powerful and conventional, although there are some advantages with using the "central" procedure.  Mainly, UCM test cannot order tables (i.e., use McNemar's conditional p-value) for a two-sided alternative with non-zero delta.  Thus, to calculate a two-sided confidence interval with UCM test, one has to resort to using the "central" approach.  For other tests, there is an equivalent statistic based on delta, and the two-sided p-value is determined by either the Agresti-Min interval (\code{tsmethod = "square"}) or Chan-Zhang interval (\code{tsmethod = "central"}).

The CSM test is computationally intensive due to iteratively maximizing the p-value calculation to order the tables.  The CSM ordering matrix has been stored for total sample sizes less than or equal to 200.  Thus, using the \code{useStoredCSM = TRUE} can greatly improve computation time.  However, the stored ordering matrix was computed with \code{npNumbers=100} and it is possible that the ordering matrix was not optimal for larger \code{npNumbers}.  Increasing \code{npNumbers} and setting \code{useStoredCSM = FALSE} ensures the p-value is correctly calculated at the expense of significantly greater computation time.  The stored ordering matrix is not used in the calculation of confidence intervals or non-inferiority tests, so CSM can still be very computationally intensive.
}
\value{
A list with class "htest" containing the following components:

\item{p.value}{The p-value of the test}
\item{statistic}{The observed test statistic to determine more extreme tables}
\item{estimate}{An estimate of the difference in proportion}
\item{null.value}{The difference in proportion under the null}
\item{conf.int}{A confidence interval for the difference in proportion.  Only present if \code{conf.int = TRUE}}
\item{alternative}{A character string describing the alternative hypothesis}
\item{method}{A character string describing the method to determine more extreme tables}
\item{tsmethod}{A character string describing the method to implement two-sided tests}
\item{np}{The nuisance parameter that maximizes the p-value}
\item{np.range}{The range of nuisance parameters considered}
\item{data.name}{A character string giving the names of the data}
\item{parameter}{The sample sizes}
}

\references{
Berger, R.L. and Sidik, K. (2003) Exact unconditional tests for 2 x 2 matched-pairs design. \emph{Statistical Methods in Medical Research}, \bold{12}, 91--108

Hsueh, H., Liu, J., and Chen, J.J. (2001) Unconditional exact tests for equivalence or noninferiority for paired binary endpoints. \emph{Biometrics}, \bold{57}, 478--483
}

\author{Peter Calhoun}
\note{
CSM test is much more computationally intensive.  Increasing the number of nuisance parameters considered and refining the p-value will increase the computation time, but more likely to ensure accurate calculations. Performing confidence intervals also greatly increases computation time.
}
\section{Warning}{CSM tests with confidence intervals may take a very long time, even for small sample sizes.}

\seealso{
\code{mcnemar.test} and \pkg{exact2x2}
}

\examples{
data <- matrix(c(3,6,1,0), 2, 2,
               dimnames=list(c("Population 1 Success", "Population 1 Failure"),
                             c("Population 2 Success", "Population 2 Failure")))
paired.exact.test(data, method="UAM", alternative="less")

\dontrun{
# Ensure that the ExactData R package is available before running the CSM test.
if (requireNamespace("ExactData", quietly = TRUE)) {
paired.exact.test(data, method="CSM", alternative="less")
}
}
}

\keyword{htest}
\concept{Barnard}
\concept{McNemar}
\concept{Unconditional}
\concept{Exact Tests}
