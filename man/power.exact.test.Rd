\name{power.exact.test}
\alias{power.exact.test}
\title{
Power calculation for unconditional exact test
}
\description{
Calculates the power of the design for known sample sizes and true probabilities.
}
\usage{
power.exact.test(p1, p2, n1, n2, npNumbers = 100, alpha = 0.05, 
                 alternative = "two.sided", interval = FALSE, beta = 0.001, 
                 method = "Z-pooled", ref.pvalue = TRUE, simulation = FALSE, 
                 nsim = 100)
}
\arguments{
  \item{p1}{
The probability of success given in group 1
}
  \item{p2}{
The probability of success given in group 2
}
  \item{n1}{
The sample size in group 1
}
  \item{n2}{
The sample size in group 2
}
  \item{npNumbers}{
Number: The number of nuisance parameters considered
}
  \item{alpha}{
Significance level
}
  \item{alternative}{
Indicates the alternative hypothesis: must be either "less", "two.sided", or "greater"
}
  \item{interval}{
Logical: Indicates if a confidence interval on the nuisance parameter should be computed
}
  \item{beta}{
Number: Confidence level for constructing the interval of nuisance parameters considered.
Only used if interval=TRUE
}
  \item{method}{
Indicates the method for finding tables as or more extreme than the observed table:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", "Boschloo", "CSM", "CSM modified", or "CSM approximate"
}
  \item{ref.pvalue}{
Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected
}
  \item{simulation}{
Logical: Indicates if the power calculation is exact or estimated by simulation
}
  \item{nsim}{
Number of simulations run.  Only used if simulation=TRUE
}
}

\details{
The power calculations are for binomial models.  The design must know the fixed sample sizes
in advance.  There are \eqn{(n_1+1) \times (n_2+1)} possible tables that could be produced.  There are two ways to
calculate the power: simulate the tables under two independent binomial distributions or consider all possible
tables and calculate the exact power.  The calculations can be done for any \code{exact.test}
computation or using Fisher's exact test.
}
\value{The function returns the computed power.}
\references{
Berger, R. (1994) Power comparison of exact unconditional tests for comparing two binomial proportions. 
Institute of Statistics Mimeo Series No. 2266

Berger, R. (1996) More powerful tests from confidence interval p values. American Statistician, 50, 314-318

Boschloo, R. D. (1970), Raised Conditional Level of Significance for the 2x2-table
when Testing the Equality of Two Probabilities. Statistica Neerlandica, 24, 1-35
}
\author{
Peter Calhoun
}
\note{
The code takes a very long time for the CSM test.  Not refining the p-value often yields similar
results and decreases the computation time.
}

\seealso{
\pkg{statmod}
}
\examples{

power.exact.test(0.20,0.80,10,20)
power.exact.test(0.20,0.80,10,20,method="Fisher")
set.seed(218461)
power.exact.test(0.20,0.80,10,20,interval=TRUE,method="Boschloo",simulation=TRUE,nsim=100)

}

\keyword{Power}
\keyword{Barnard}
\keyword{Boschloo}
\keyword{Unconditional}
\keyword{Exact}
