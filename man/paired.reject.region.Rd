\name{paired.reject.region}
\alias{paired.reject.region}
\title{Rejection Region for 2x2 Tables with Paired Samples}
\description{Determines the rejection region for comparing two paired proportions.}
\usage{
paired.reject.region(N, alternative = c("two.sided", "less", "greater"),
    alpha = 0.05, npNumbers = 100, np.interval = FALSE, beta = 0.001,
    method = c("uam", "ucm", "uamcc", "csm", "cm", "am", "amcc"),
    tsmethod = c("square", "central"),
    delta = 0, convexity = TRUE, useStoredCSM = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{The total sample size}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{alpha}{Significance level}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding the more extreme tables:
must be either "UAM", "UCM", "UAMCC", "CSM", "CM", "AM", or "AMCC"}
  \item{tsmethod}{A character string describing the method to implement two-sided tests}
  \item{delta}{Number: null hypothesis of the difference in proportion}
  \item{convexity}{Logical: assumes convexity for interval approach.  Only used if np.interval=TRUE}
  \item{useStoredCSM}{Logical: uses stored CSM ordering matrix.  Only used if method="csm"}
}
\details{
The rejection region is calculated for paired samples.  Rejection region can be determined for any unconditional exact test in \code{\link{paired.exact.test}}, the Conditional McNemar's (CM) exact test, the Asymptotic McNemar's (AM) test, or Asymptotic McNemar's test with Continuity Correction (AMCC) (note: asymptotic tests are not exact tests).  In very rare cases, using the nuisance parameter interval approach does not attain the convexity property, so it is possible using \code{convexity=TRUE} could yield an inaccurate power calculation with this method.  This is \emph{extremely unlikely} though, so default is to assume convexity and speed up computation time.  For details regarding parameters, see \code{\link{paired.exact.test}}.
}
\value{
A matrix of the rejection region.  The rows and columns represent the discordant pairs.  Specifically, the columns represent x12, the number of successes in first group and number of failures in second group, and rows represent x21, the number of failures in first group and number of successes in second group.  The number of concordant pairs is simply the total sample size minus number of discordant pairs.  The cells represent whether the test is rejected (1) or failed to be rejected (0).  Values with an NA are not possible.  This matrix represents all possible 2x2 tables.
}
\author{Peter Calhoun}
\note{
McNemar's asymptotic tests are not exact test and may have inflated type 1 error rates.  These options were added to compute the rejection region efficiently when using asymptotic tests.
}
\seealso{\code{\link{power.paired.test}}}
\examples{
\dontrun{
# Ensure that the ExactData R package is available before running the CSM test.
if (requireNamespace("ExactData", quietly = TRUE)) {
paired.reject.region(N=15, alternative="two.sided", method="CSM")
}
}

paired.reject.region(N=15, alternative="less", method="UAM", delta=0.10)
}
\keyword{htest}
\concept{Rejection region}
\concept{Barnard}
\concept{McNemar}
\concept{Unconditional}
\concept{Exact Tests}
