\name{p.value.compute}
\alias{p.value.compute}
\title{Computation of the p-value for a given hypotheses test}
\description{
 Computes the p-values for all the entries in the matrix \code{test.mat} using the asymtotic properties of the estimator of the wavelet correlation given in (Whitcher, 2000). 
}
\usage{
p.value.compute(test.mat, var.ind.mat = 0, n.ind = 0, test.method = "gaussian",
                 proc.length, sup, num.levels, use.tanh = FALSE)
}
\arguments{
  \item{test.mat}{ matrix containing the wavelet correlation to be tested}
\item{var.ind.mat}{ matrix containing the variance inter individuals of the correlation. Only used with \code{test.method="t.test"}. (default not used)}
  \item{n.ind}{number of individuals to take into account in the test. Only used with \code{test.method="t.test"}. (default not used)}  
 \item{test.method}{ name of the method to be applied. \code{"gaussian"} assumes a gaussian law for the estimator. \code{"t.test"} implements a t.test for computing the p-value. (default \code{"gaussian"}) }
  \item{proc.length}{ specifies the length of the original processes using to construct the \code{cor.mat}}
\item{num.levels}{ specifies the number of the wavelet scale to take into account in the hypothesis test. Only used with \code{test.method="gaussian"} }
  \item{use.tanh}{ logical. If FALSE take the \code{atanh} of the correlation values before applying the hypothesis test, in order to use the Fisher approximation}  
   \item{sup}{ indicates the correlation threshold to consider in each hypothesis test.
}
}
 \details{
 Each hypothesis test is written as :
\eqn{H_0} : "|correlation| \eqn{<=} sup"
\eqn{H_1} : "|correlation| \eqn{>} sup"
This function is essentially an internal function called by \code{const.adj.mat}.
}

\value{
 Vector with the p-value for each entry of the matrix.
}
\references{
 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
}
\author{S. Achard}
\seealso{code{\link{const.adj.mat}}}
\examples{
data(brain)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]


# Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 4, 
                               boundary = "periodic", p.corr = 0.975)

# For scale 4
pvalue.cor<-p.value.compute(wave.cor.list[[4]],proc.length=dim(brain)[1], sup=0.44, 
                            num.levels=4)

}
\keyword{multivariate}
\keyword{htest}
