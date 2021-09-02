\name{exact.reject.region}
\alias{exact.reject.region}
\title{Rejection Region for 2x2 Tables with Independent Samples}
\description{Determines the rejection region for comparing two independent proportions.}

\usage{
exact.reject.region(n1, n2, alternative = c("two.sided", "less", "greater"),
    alpha = 0.05, npNumbers = 100, np.interval = FALSE, beta = 0.001,
    method = c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
               "csm", "fisher", "pearson chisq", "yates chisq"),
    tsmethod = c("square", "central"), delta = 0, convexity = TRUE,
    useStoredCSM = TRUE)
}

\arguments{
  \item{n1}{The sample size in first group}
  \item{n2}{The sample size in second group}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{alpha}{Significance level}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding the more extreme tables:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", "Boschloo", "CSM", "Fisher", "Pearson Chisq", or "Yates Chisq"}
  \item{tsmethod}{Indicates two-sided method: must be either "square" or "central"}
  \item{delta}{Number: null hypothesis of the difference in proportion}
  \item{convexity}{Logical: assumes convexity for interval approach.  Only used if np.interval=TRUE}
  \item{useStoredCSM}{Logical: uses stored CSM ordering matrix.  Only used if method="csm"}
}

\details{
The rejection region is calculated for binomial models with independent samples.  The design must know the fixed sample sizes in advance.  Rejection region can be determined for any unconditional exact test in \code{\link{exact.test}}, Fisher's exact test, or chi-square test (Yates' or Pearson's; note: these are not exact tests).  In very rare cases, using the nuisance parameter interval approach does not attain the convexity property, so it is possible using \code{convexity=TRUE} could yield an inaccurate power calculation with this method.  This is \emph{extremely unlikely} though, so default is to assume convexity and speed up computation time.  For details regarding parameters, see \code{\link{exact.test}}.
}
\value{
A matrix of the rejection region.  The columns represent the number of successes in first group, rows represent the number of successess in second group, and cells represent whether the test is rejected (1) or failed to be rejected (0).  This matrix represents all possible 2x2 tables.
}
\references{
Barnard, G.A. (1947) Significance tests for 2x2 tables. \emph{Biometrika}, \bold{34}, 123-138

Chan, I. (2003), Proving non-inferiority or equivalence of two treatments with dichotomous endpoints using exact methods, \emph{Statistical Methods in Medical Research}, \bold{12}, 37--58
}
\author{Peter Calhoun}
\note{
Pearson's and Yates' chi-square tests are not exact tests, so the function name may be a misnomer.  These tests may have inflated type 1 error rates.  These options were added to compute the rejection region efficiently when using asymptotic tests.
}
\seealso{\code{\link{power.exact.test}}}
\examples{

\dontrun{
# Ensure that the ExactData R package is available before running the CSM test.
if (requireNamespace("ExactData", quietly = TRUE)) {
exact.reject.region(n1=10, n2=20, alternative="two.sided", method="CSM")
}
}

exact.reject.region(n1=10, n2=20, alternative="less", method="Z-pooled", delta=0.10)
}
\keyword{htest}
\concept{Rejection region}
\concept{Barnard}
\concept{Boschloo}
\concept{Unconditional}
\concept{Exact Tests}
