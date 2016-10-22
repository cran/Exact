\name{exact.test}
\alias{exact.test}
\alias{binomialCode}
\alias{multinomialCode}
\title{Unconditional exact tests for 2x2 tables}
\description{Calculates Barnard's or Boschloo's unconditional exact test for binomial or multinomial models}
\usage{
exact.test(data, alternative = c("two.sided", "less", "greater"), npNumbers = 100,
           beta = 0.001, interval = FALSE, 
             method = c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
                      "csm", "csm approximate", "csm modified"),
              model = c("Binomial", "Multinomial"), cond.row = TRUE, to.plot = TRUE,
           ref.pvalue = TRUE)
}

\arguments{
  \item{data}{A two dimensional contingency table in matrix form}
  \item{alternative}{Indicates the alternative hypothesis: must be either "less", "two.sided", or "greater"}
  \item{npNumbers}{Number: The number of nuisance parameters considered}
  \item{beta}{Number: Confidence level for constructing the interval of nuisance parameters considered.  Only used if interval=TRUE}
  \item{interval}{Logical: Indicates if a confidence interval on the nuisance parameter should be computed}
  \item{method}{
Indicates the method for finding tables as or more extreme than the observed table:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", "Boschloo", "CSM", "CSM approximate", or "CSM modified".
CSM tests cannot be calculated for multinomial models
}
  \item{model}{The model being used: must be either "Binomial" or "Multinomial"}
  \item{cond.row}{Logical: Indicates if row margins are fixed in the binomial models.  Only used if model="Binomial"}
  \item{to.plot}{Logical: Indicates if plot of p-value vs. nuisance parameter should be generated.  Only used if model="Binomial"}
  \item{ref.pvalue}{Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected.  Only used if model="Binomial"}
}

\details{
Unconditional exact tests can be used for binomial or multinomial models.  The binomial model assumes the row or column margins (but not both) are known in advance, while the multinomial model assumes only the total sample size is known beforehand.  Conditional tests have both row and column margins fixed.  The null hypothesis is that the rows and columns are independent. Under the binomial model, the user will need to input which margin is fixed (default is rows).

See the following formulas in the Referene Manual: \url{https://CRAN.R-project.org/package=Exact}.

Let \eqn{X} denote a generic 2x2 table with fixed sample sizes \eqn{n_1} and \eqn{n_2}, \eqn{X_0} denote the observed table, and \eqn{T(X)} represent the test statistic function.  The null hypothesis can be written as \eqn{p_1=p_2 \equiv p}. The p-value function with rows fixed is the product of two independent binomials:

\deqn{P(X|p)= \sup_{0 \leq p \leq 1} \sum_{T(X) \geq T(X_0)} {n_1 \choose x_{11}} {n_2 \choose x_{21}} p^{x_{11}+x_{21}} (1-p)^{x_{12}+x_{13}}}

The multinomial model is similar except the summand has a multinomial distribution with two nuisance parameters.

There are several possible test statistics to determine the 'as or more extreme' tables seen in the index of summation.  The \code{method} variable lets the user choose the test statistic being used.  A brief description for each test statistic is given below (see References for more details):

Let \eqn{\hat{p_1}=x_{11}/n_1}, \eqn{\hat{p_2}=x_{21}/n_2}, and \eqn{\hat{p}=(x_{11}+x_{21})/(n_1+n_2)}.

Z-unpooled (or Wald):
\deqn{Z_u(x_{11},x_{21})=\frac{\hat{p_2}-\hat{p_1}}{\sqrt{\frac{\hat{p_1}(1-\hat{p_1})}{n_1}+\frac{\hat{p_2}(1-\hat{p_2})}{n_2}}}}

Z-pooled (or Score):
\deqn{Z_p(x_{11},x_{21})=\frac{\hat{p_2}-\hat{p_1}}{\sqrt{\frac{\hat{p}(1-\hat{p})}{n_1}+\frac{\hat{p}(1-\hat{p})}{n_2}}}}

Santner and Snell:
\deqn{D(x_{11},x_{21})=\hat{p_2}-\hat{p_1}}

Boschloo:

Uses the p-value from Fisher's exact test as the test statistic.

CSM:

Starts with the most extreme table and adds other 'as or more extreme' tables one step at a time by maximizing the summand of the p-value function.  This approach can be computationally intensive.

CSM modified:

Starts with all tables that must be more extreme and adds other 'as or more extreme' tables one step at a time by maximizing the summand of the p-value function.  This approach can be computationally intensive.

CSM approximate:

Maximizes the summand of the p-value function for each possible table.  Thus, the test statistic is the p-value function without the summation.  This approach is less computationally intensive than the CSM test because the maximization is not repeated at each step. 

The supremum of the common success probability is taken over all values between 0 and 1.  Another approach, proposed by Berger and Boos, is to take the supremum over a Clopper-Pearson confidence interval.  This approach adds a small penalty to the p-value to ensure a level-\eqn{\alpha} test, but eliminates unlikely probabilities from inflating the p-value.  The p-value function becomes:

\deqn{P(X|p)= \left(\sup_{p \in C_\beta} \sum_{T(X) \geq T(X_0)} {n_1 \choose x_{11}} {n_2 \choose x_{21}} p^{x_{11}+x_{21}} (1-p)^{x_{12}+x_{13}}\right) + \beta }

where \eqn{C_\beta} is the \eqn{100(1-\beta)\%} confidence interval of \eqn{p}

There are many ways to define the two-sided p-value; this code uses the \code{fisher.test} approach by summing the probabilities for both sides of the table.
}

\value{
A list with class "htest" containing the following components:

\item{p.value}{The computed p-value}
\item{test.statistic}{The observed test statistic}
\item{estimate}{An estimate of the parameter tested}
\item{alternative}{A character string describing the alternative hypothesis}
\item{model}{A character string describing the model design ("Binomial" or "Multinomial")}
\item{method}{A character string describing the method to determine 'as or more extreme' tables}
\item{np}{The nuisance parameter that maximizes the p-value.  For multinomial models, both nuisance parameters are given}
\item{np.range}{The range of nuisance parameters considered.  For multinomial models, both nuisance parameter ranges are given}
\item{data.name}{A character string giving the names of the data}
}
\references{This code was influenced by the FORTRAN program located at \url{http://www4.stat.ncsu.edu/~boos/exact/}}
\author{Peter Calhoun}
\note{
CSM test and multinomial models are much more computationally intensive.  I have also spent a greater amount of time making the computations for the binomial models more efficient; future work will be devoted to improving the multinomial models. Boschloo's test also takes longer due to calculating Fisher's p-value for every possible table; however, a created function that calculates Fisher's test efficiently is utilized. Increasing the number of nuisance parameters considered and refining the p-value will increase the computation time.
}

\section{Warning}{Multinomial models and CSM tests may take a very long time, even for sample sizes less than 100.}
\seealso{
\code{fisher.test} and \pkg{exact2x2}
}
\examples{
data <- matrix(c(7, 8, 12, 3), 2, 2, byrow=TRUE)
exact.test(data, alternative="less",to.plot=TRUE)
exact.test(data, alternative="two.sided", interval=TRUE, beta=0.001, npNumbers=100,
           method="Z-pooled",to.plot=FALSE)
exact.test(data, alternative="two.sided", interval=TRUE, beta=0.001, npNumbers=100,
           method="Boschloo", to.plot=FALSE)

#Example from Barnard's (1947) appendix:
data <- matrix(c(4, 0, 3, 7), 2, 2,
               dimnames=list(c("Box 1","Box 2"), c("Defective","Not Defective")))
exact.test(data, method="CSM", alternative="two.sided")

data <- matrix(c(6, 8, 4, 3), 2, 2, byrow=TRUE)
exact.test(data, model="Multinomial", alternative="less", method="Z-pooled")
}
\keyword{Barnard}
\keyword{Boschloo}
\keyword{Unconditional}
\keyword{Exact Tests}
