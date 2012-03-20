\name{exact.test}
\alias{exact.test}
\alias{Barnard}
\alias{exact2x2}
\title{Unconditional exact tests for 2x2 tables}
\description{
Calculates Barnard's or Boschloo's unconditional exact test for binomial or multinomial models
}
\usage{
exact.test(data, alternative = "less", npNumbers = 0.01, beta = 0.001, interval = FALSE,
           method = "Z-pooled", model = "Binomial")
}
\arguments{
  \item{data}{
A two dimensional contingency table in matrix form
}
  \item{alternative}{
Indicates the alternative hypothesis: must be either "less", "two.sided", or "greater"
}
  \item{npNumbers}{
Number: Increment for the search in the nuisance parameter
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
must be either "Z-pooled", "Z-unpooled", or "Boschloo"
}
  \item{model}{
The model being used: must be either "Binomial" or "Multinomial"
}
}
\details{
Unconditional exact tests can be used for binomial or multinomial models.  The binomial model assumes the rows 
or columns are known in advance, while the multinomial model only assumes the total sample size is known beforehand.  
Conditional tests have both rows and columns fixed.  Fisher's test conditions on both margins to avoid estimating a 
nuisance parameter.  Barnard's test considers all possible values for the nuisance parameter and chooses the one 
that maximizes the p-value.

There are several statistics used to define tables as more extreme than the observed table.
The method variable lets the user choose the test statistic being used.

The null hypothesis is that the rows and columns are independent.  The function assumes \emph{the rows are the fixed margins for the binomial model}.
The null hypothesis can be interpreted as the probability of the top-left cell given
the first row sample size is equal to the probability of the bottom-left cell given the second row sample size.
There are many ways to define the two.sided p-value; this version uses the fisher.test approach by summing the 
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
This code is similar to the FORTRAN program located at \url{www.stat.ncsu.edu/exact}.
}
\author{
Peter Calhoun
}
\note{
I am indebted to Roger Berger for helping me with some of the theory, Kamil Erguler for changing my one-sided
p-value calculation, Tal Galili for making my first version of the code accessible to the R community, and
Philo Calhoun for initiating and supporting my research.
}

\section{Warning}{Binomial models tend to have memory issues for sample sizes greater than 1000.
Multinomial models are much more computational intensive due to two nuisance parameters being maximized.
Boschloo's test also takes longer due to calculating Fisher's p-value for every possible table.
I suggest keeping the sample sizes small unless the user does not mind waiting.
}

\seealso{
\code{\link{fisher.test},\link{Barnard},\link{exact2x2}}
}


\examples{
data<-matrix(c(7,12,8,3),2,2)
exact.test(data,alternative="less")
exact.test(data,alternative="two.sided",interval=TRUE,beta=.001,npNumbers=0.01,method="Z-pooled")
exact.test(data,alternative="two.sided",interval=TRUE,beta=.001,npNumbers=0.01,method="Boschloo")

data<-matrix(c(6,4,8,3),2,2)
exact.test(data,model="Multinomial",alternative="less",method="Z-pooled")
}

\keyword{Barnard}
\keyword{Boschloo}
\keyword{Unconditional}
\keyword{Nonparametric}

