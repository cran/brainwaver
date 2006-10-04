\name{choose.thresh.nbedges}
\alias{choose.thresh.nbedges}
\title{Threshold associated to a given number of edges.}
\description{
Computes the threshold for the correlation matrix in order to obtain an adjacency matrix with a given number of edges.
}
\usage{
choose.thresh.nbedges(cor.mat,  var.ind.mat = 0, n.ind = 0,  thresh = 0.05,  
                       nb.edges = 405, test.method = "gaussian",  
                       proc.length = 518, num.levels, use.tanh = FALSE, 
                       max.iter = 10)
}
\arguments{
  \item{cor.mat}{ matrix containing the correlation values. (must be diagonal with 1 on the diagonal)}
  \item{var.ind.mat}{matrix containing the variance inter individuals of the correlation. Only used with \code{test.method="t.test"}. (default not used)}
  \item{n.ind}{number of individuals to take into account in the test. Only used with \code{test.method="t.test"}. (default not used)}
   \item{thresh}{ indicates the rate at which the FDR procedure is controlled. (default 0.05)}
  \item{nb.edges}{ indicates the exact number of edges that the final graph should contain.
}
  \item{test.method}{ name of the method to be applied. \code{"gaussian"} assumes a gaussian law for the estimator. \code{"t.test"} implements a t.test for computing the p-value. (default \code{"gaussian"}) }
  \item{proc.length}{ specifies the length of the original processes using to construct the \code{cor.mat}}
\item{num.levels}{ specifies the number of the wavelet scale to take into account in the hypothesis test. Only used with \code{test.method="gaussian"} }
  \item{use.tanh}{ logical. If FALSE take the \code{atanh} of the correlation values before applying the hypothesis test, in order to use the Fisher approximation}
 \item{max.iter}{indicates the number of maximum iteration to compute before 
stopping the loop}
}
\details{
In order to compare graphs, the best way to do it is to make sure that all the graphs have the same number of edges!  
}
\value{
Real number corresponding to the threshold value.
}
\references{
 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
}
\author{S. Achard}
\note{only in version 2 and higher}
\seealso{}
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
nb.edges<-sum(adj.mat.4)/2

sup.thresh<-choose.thresh.nbedges(wave.cor.list[[4]],nb.edges=nb.edges,
                         proc.length=dim(brain)[1],num.levels=4)

}
\keyword{graphs}
