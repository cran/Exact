\name{power.exact.test}
\alias{power.exact.test}
\title{Power Calculations for 2x2 Tables}
\description{Calculates the power of the design for known sample sizes and true probabilities.}
\usage{
power.exact.test(p1, p2, n1, n2, alternative = c("two.sided", "less", "greater"),
    alpha = 0.05, npNumbers = 100, np.interval = FALSE, beta = 0.001,
    method = c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
               "csm", "csm approximate", "fisher", "chisq", "yates chisq"),
    tsmethod = c("square", "central"), ref.pvalue = TRUE, simulation = FALSE,
    nsim = 100, delta = 0, convexity = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{The probability of success given in first group}
  \item{p2}{The probability of success given in second group}
  \item{n1}{The sample size in first group}
  \item{n2}{The sample size in second group}
  \item{alternative}{Indicates the alternative hypothesis: must be either "two.sided", "less", or "greater"}
  \item{alpha}{Significance level}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{np.interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if np.interval=TRUE}
  \item{method}{
Indicates the method for finding more extreme tables:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", "Boschloo", "CSM", "CSM approximate", "Fisher", "Chisq", or "Yates Chisq"
}
  \item{tsmethod}{Indicates two-sided method: must be either "square" or "central"}
  \item{ref.pvalue}{Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected}
  \item{simulation}{Logical: Indicates if the power calculation is exact or estimated by simulation}
  \item{nsim}{Number of simulations run.  Only used if simulation=TRUE}
  \item{delta}{Number: null hypothesis of the difference in proportion.  Only used if method is "z-pooled" or "csm"}
  \item{convexity}{Logical: assumes convexity for interval approach.  Only used if np.interval=TRUE}
}
\details{
The power calculations are for binomial models.  The design must know the fixed sample sizes in advance.  There are (n1+1) x (n2+1) possible tables that could be produced.  There are two ways to calculate the power: simulate the tables under two independent binomial distributions or determine the rejection region for all possible tables and calculate the exact power.  The calculations can be done using any \code{\link{exact.test}} computation, Fisher's exact test, or chi-square tests (Yates' or Pearson's; note: these are not exact tests).  The power calculations utilize the convexity property, which greatly speeds up computation time (see \code{\link{exact.reject.region}} documentation).
}
\value{
A list with class "power.htest" containing the following components:

\item{n1, n2}{The respective sample sizes}
\item{p1, p2}{The respective proportions}
\item{alpha}{Significance level}
\item{power}{Power of the test}
\item{alternative}{A character string describing the alternative hypothesis}
\item{delta}{Null hypothesis of the difference in proportion}
\item{method}{A character string describing the method to determine more extreme tables}
}
\references{
Berger, R. (1994) Power comparison of exact unconditional tests for comparing two binomial proportions. \emph{Institute of Statistics Mimeo Series No. 2266}

Berger, R. (1996) More powerful tests from confidence interval p values. \emph{American Statistician}, \bold{50}, 314-318

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table when Testing the Equality of Two Probabilities. \emph{Statistica Neerlandica}, \bold{24}, 1-35
}
\author{Peter Calhoun}
\note{
Pearson's and Yates' chi-square tests are not exact test, so the function name may be a misnomer.  However, these options were added to compute the power efficiently when using asymptotic tests.
}
\seealso{\code{\link{exact.reject.region}} and \pkg{statmod}}
\examples{
# Superiority power #
power.exact.test(p1=0.15, p2=0.60, n1=15, n2=30, method="Z-pooled")
power.exact.test(p1=0.15, p2=0.60, n1=15, n2=30, method="CSM")
power.exact.test(p1=0.15, p2=0.60, n1=15, n2=30, method="Fisher")
power.exact.test(p1=0.15, p2=0.60, n1=15, n2=30, method="Boschloo",
                 np.interval=TRUE, beta=0.001)

# Non-inferiority power #
power.exact.test(p1=0.30, p2=0.30, n1=65, n2=65, method="Z-pooled",
                 delta=0.2, alternative="less")
}
\keyword{htest}
\concept{Power}
\concept{Barnard}
\concept{Boschloo}
\concept{Unconditional}
\concept{Exact Tests}
