\name{exact.test}
\alias{exact.test}
\title{Unconditional exact tests for 2x2 tables}
\description{
Calculates Barnard's or Boschloo's unconditional exact test for binomial or multinomial models
}
\usage{
exact.test(data, alternative = "two.sided", npNumbers = 100, beta = 0.001,
           interval = FALSE, method = "Z-pooled", model = "Binomial", 
           cond.row = TRUE, to.plot = TRUE, ref.pvalue=TRUE)
}
\arguments{
  \item{data}{
A two dimensional contingency table in matrix form
}
  \item{alternative}{
Indicates the alternative hypothesis: must be either "less", "two.sided", or "greater"
}
  \item{npNumbers}{
Number: The number of nuisance parameters considered
}
  \item{beta}{
Number: Confidence level for constructing the interval of nuisance parameters considered.
Only used if interval=TRUE
}
  \item{interval}{
Logical: Indicates if a confidence interval on the nuisance parameter should be computed
}
  \item{method}{
Indicates the method for finding tables as or more extreme than the observed table:
must be either "Z-pooled", "Z-unpooled", "Santner and Snell", or "Boschloo"
}
  \item{model}{
The model being used: must be either "Binomial" or "Multinomial"
}
  \item{cond.row}{
Logical: Indicates if row margins are fixed in the binomial models.  Only used if model="Binomial"
}
  \item{to.plot}{
Logical: Indicates if plot of p-value vs nuisance parameter should be generated.  Only used if model="Binomial"
}
  \item{ref.pvalue}{
Logical: Indicates if p-value should be refined by maximizing the p-value function after the nuisance parameter is selected.  Only used if model="Binomial"
}
}
\details{
Unconditional exact tests can be used for binomial or multinomial models.  The binomial model assumes the row 
or column margins (but not both) are known in advance, while the multinomial model only assumes the total sample size is known beforehand.  
Conditional tests have both row and column margins fixed.  Fisher's test conditions on both margins to avoid estimating a 
nuisance parameter.  Barnard's test considers all possible values for the nuisance parameter and chooses the one 
that maximizes the p-value.

There are several statistics used to define tables as or more extreme than the observed table.
The method variable lets the user choose the test statistic being used.

The null hypothesis is that the rows and columns are independent.  Under the binomial model, the user will need to input which margin is 
fixed (default is rows).  The null hypothesis can be interpreted as the probability of success given in the first group is equal to the probability 
of success given in the second group. There are many ways to define the two.sided p-value; this code uses the fisher.test() approach by summing the 
probabilities for both sides of the table.
}
\value{
\item{p.value}{The computed p-value}
\item{test.statistic}{The observed test statistic}
\item{np}{The nuisance parameter that maximizes the p-value.  For multinomial models, both nuisance parameters
are given}
\item{np.range}{The range of nuisance parameters considered.  For multinomial models, both nuisance parameter
ranges are given}
}
\references{
This code was influenced by the FORTRAN program located at \url{www.stat.ncsu.edu/exact}
}
\author{
Peter Calhoun
}
\note{
Multinomial models are much more computationally intensive due to more tables considered and maximizing over two 
nuisance parameters.  I also have spent a greater amount of time making the computations for the binomial models more efficient;
future work will be devoted to improving the multinomial models.  Boschloo's test also takes longer due to
calculating Fisher's p-value for every possible table; however, a created function that calculates Fisher's test efficiently
is utilized. Increasing the number of nuisance parameters considered and refining the p-value will increase the computation time.
}

\section{Warning}{Multinomial models may take a very long time, even for sample sizes less than 100.}

\seealso{
\code{\link{fisher.test}}, \code{Barnard}, and \code{exact2x2}
}
\examples{
data<-matrix(c(7,8,12,3),2,2,byrow=TRUE)
exact.test(data,alternative="less",to.plot=TRUE)
exact.test(data,alternative="two.sided",interval=TRUE,beta=0.001,npNumbers=100,method="Z-pooled",
           to.plot=FALSE)
exact.test(data,alternative="two.sided",interval=TRUE,beta=0.001,npNumbers=100,method="Boschloo",
           to.plot=FALSE)

data<-matrix(c(6,8,4,3),2,2,byrow=TRUE)
exact.test(data,model="Multinomial",alternative="less",method="Z-pooled")
}

\keyword{Barnard}
\keyword{Boschloo}
\keyword{Unconditional}
\keyword{Exact}
\keyword{Nonparametric}

