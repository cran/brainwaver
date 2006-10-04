\name{const.adj.mat}
\alias{const.adj.mat}
\title{Computation of the adjacency matrix}
\description{
 Computes the adjacency matrix for a given correlation matrix.
}
\usage{
const.adj.mat(cor.mat, var.ind.mat = 0, n.ind = 0, thresh = 0.05, sup = 0, 
              test.method = "gaussian", proc.length, num.levels, 
              use.tanh = FALSE)
}

\arguments{
  \item{cor.mat}{ matrix containing the correlation values. (must be diagonal with 1 on the diagonal)}
  \item{var.ind.mat}{matrix containing the variance inter individuals of the correlation. Only used with \code{test.method="t.test"}. (default not used)}
  \item{n.ind}{number of individuals to take into account in the test. Only used with \code{test.method="t.test"}. (default not used)}
   \item{thresh}{ indicates the rate at which the FDR procedure is controlled. (default 0.05)}
  \item{sup}{ indicates the correlation threshold to consider in each hypothesis test.
}
  \item{test.method}{ name of the method to be applied. \code{"gaussian"} assumes a gaussian law for the estimator. \code{"t.test"} implements a t.test for computing the p-value. (default \code{"gaussian"}) }
  \item{proc.length}{ specifies the length of the original processes using to construct the \code{cor.mat}}
\item{num.levels}{ specifies the number of the wavelet scale to take into account in the hypothesis test. Only used with \code{test.method="gaussian"} }
  \item{use.tanh}{ logical. If FALSE take the \code{atanh} of the correlation values before applying the hypothesis test, in order to use the Fisher approximation}
}
 \details{
 Each hypothesis test is written as :

\eqn{H_0} : "|correlation| \eqn{<=} sup"

\eqn{H_1} : "|correlation| \eqn{>} sup"

}
\value{
 Binary matrix.
}

\references{
 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
}
\author{S. Achard}
\note{ }
\seealso{ \code{\link{const.adj.list}}}
\examples{
data(brain)
brain<-as.matrix(brain)
# WARNING : To process only the first five regions
brain<-brain[,1:5]


#Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 6,
                               boundary = "periodic", p.corr = 0.975)

#Construction of the adjacency matrice for scale 4

adj.mat.4<-const.adj.mat(wave.cor.list[[4]], sup = 0.44,proc.length=dim(brain)[1], 
                          num.levels=4)

image(adj.mat.4,col=gray((0:20)/20))
}
\keyword{graphs}
\keyword{ts}
\keyword{multivariate}

