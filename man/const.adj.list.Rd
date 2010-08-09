\name{const.adj.list}
\alias{const.adj.list}
\title{Computation of the list of adjacency matrices}
\description{
Computes the list of the adjacency matrices in terms of the scale of the wavelet decomposition.
}
\usage{
const.adj.list(wave.cor.list, wave.var.ind = 0, n.ind = 0, thresh = 0.05, 
                 sup = 0, test.method = "gaussian", proc.length, 
                 use.tanh = FALSE)
}

\arguments{
  \item{wave.cor.list}{object of class \code{"Wave Correlation"} containing the correlation matrices to be analysed}
  \item{wave.var.ind}{object of class \code{"Wave Correlation"} containing the inter individuals variance of the correlation. Only used with \code{test.method="t.test"}. (default not used)}
  \item{n.ind}{ number of individuals to take into account in the test. Only used with \code{test.method="t.test"}. (default not used)}
  \item{thresh}{ indicates the rate at which the FDR procedure is controlled. (default 0.05)}
  \item{sup}{ indicates the correlation threshold to consider in each hypothesis test
}
  \item{test.method}{ name of the method to be applied. \code{"gaussian"} assumes a gaussian law for the estimator. \code{"t.test"} implements a t.test for computing the p-value. (default \code{"gaussian"}) .}
  \item{proc.length}{ specifies the length of the original processes using to construct the \code{wave.cor.list}}
  \item{use.tanh}{ logical. If FALSE take the \code{atanh} of the correlation values before applying the hypothesis test, in order to use the Fisher approximation}
}
\details{
Each hypothesis test is written as :

\eqn{H_0} : "|correlation| \eqn{<=} sup"

\eqn{H_1} : "|correlation| \eqn{>} sup"

}
\value{
 Object of class \code{"Wave Adjacency matrix"}, basically, a list with the following
  components
 \item{d?}{Adjacency matrix for each scale of the wavelet decomposition}

 }
\references{
 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
}
\author{S. Achard}
\seealso{ \code{\link{const.adj.mat}}}
\examples{
data(brain)
brain<-as.matrix(brain)
# WARNING : To process only the first five regions
brain<-brain[,1:5]


# Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 6, 
                               boundary = "periodic", p.corr = 0.975)

# Construction of the adjacency matrices associated to each level of the 
# wavelet decomposition
wave.adj.list<-const.adj.list(wave.cor.list, sup = 0.44, proc.length=dim(brain)[1])

par(mfrow=c(3,2))

for(i in 1:4)
{
name.txt<-paste("Level ",i,sep="")
image(wave.adj.list[[i]],col=gray((0:20)/20),main=name.txt)
}

}
\keyword{graphs}
\keyword{ts}
\keyword{multivariate}