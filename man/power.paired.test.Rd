\name{power.paired.test}
\alias{power.paired.test}
\title{Power Calculations for 2x2 Tables with Paired Samples}
\description{Calculates the power of the design for known sample size and true probabilities.}
\usage{
power.paired.test(p12, p21, N, alternative = c("two.sided", "less", "greater"),
    alpha = 0.05, npNumbers = 100, np.interval = FALSE, beta = 0.001,
    method = c("mcnemar", "mcnemar with cc", "csm", "conditional exact mcnemar", 
               "asymptotic mcnemar", "asymptotic mcnemar with cc"),
    simulation = FALSE, nsim = 100,
    delta = 0, convexity = TRUE, useStoredCSM = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p12}{The probability of success in first group and failure in second group.  This is the probability of the discordant pair x12}
  \item{p21}{The probability of failure in first group and success in second group.  This is the probability of the discordant pair x21}
  \item{N}{The total sample size}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{alpha}{Significance level}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding the more extreme tables:
must be either "McNemar", "McNemar with CC", "CSM", "Conditional Exact McNemar", "Asymptotic McNemar", or "Asymptotic McNemar with CC"
}
  \item{simulation}{Logical: Indicates if the power calculation is exact or estimated by simulation}
  \item{nsim}{Number of simulations run.  Only used if simulation=TRUE}
  \item{delta}{Number: null hypothesis of the difference in proportion.  Cannot be used if method="Conditional Exact McNemar"}
  \item{convexity}{Logical: assumes convexity for interval approach.  Only used if np.interval=TRUE}
  \item{useStoredCSM}{Logical: uses stored CSM ordering matrix.  Only used if method="csm"}
}
\details{
The power calculations are for paired samples.  All possible tables can be represented by an (N+1) x (N+1) matrix.  There are two ways to calculate the power: simulate the tables under a trinomial distribution or determine the rejection region for all possible tables and calculate the exact power.  The calculations can be done using any \code{\link{paired.exact.test}} computation, McNemar's conditional exact test, or McNemar's asymptotic test with or without continuity correction (note: these are not exact tests).  The power calculations utilize the convexity property, which greatly speeds up computation time (see \code{\link{paired.reject.region}} documentation).
}
\value{
A list with class "power.htest" containing the following components:

\item{N}{The total sample size}
\item{p12, p21}{The respective discordant probabilities}
\item{alpha}{Significance level}
\item{power}{Power of the test}
\item{alternative}{A character string describing the alternative hypothesis}
\item{delta}{Null hypothesis of the difference in proportion}
\item{method}{A character string describing the method to determine more extreme tables}
}
\references{
Berger, R.L. and Sidik, K. (2003) Exact unconditional tests for 2 x 2 matched-pairs design. \emph{Statistical Methods in Medical Research}, \bold{12}, 91--108
}
\author{Peter Calhoun}
\note{
McNemar's asymptotic tests are not exact test and may have inflated type 1 error rates.  These options were added to compute the power efficiently when using asymptotic tests.
}
\seealso{\code{\link{paired.reject.region}}}

\examples{

# Superiority power #
power.paired.test(p12=0.15, p21=0.45, N=40, method="McNemar")
\dontrun{
# Ensure that the ExactData R package is available before running the CSM test.
if (requireNamespace("ExactData", quietly = TRUE)) {
power.paired.test(p12=0.15, p21=0.45, N=40, method="CSM")
}
}

# Non-inferiority power #
power.paired.test(p12=0.30, p21=0.30, N=80, method="McNemar",
                  alternative="less", delta=0.2)
}
\keyword{htest}
\concept{Power}
\concept{Barnard}
\concept{McNemar}
\concept{Unconditional}
\concept{Exact Tests}
