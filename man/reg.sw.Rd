\name{reg.sw}
\alias{reg.sw}
\title{Small-world parameters for a lattice graph}
\description{
Computes the degree, the minimum path length and the clustering coefficient for a lattice graph.
}
\usage{
reg.sw(n.nodes.rand, n.edges.rand, dist.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n.nodes.rand}{number of nodes of the graph}
  \item{n.edges.rand}{number of edges of the graph}
\item{dist.mat}{matrix with a distance associated to each pair of nodes of the graph to take into account in the computation of the small-world parameters.}
}
\value{
\item{in.degree}{mean of the degree for the whole graph.}
  \item{Lp.rand}{mean of the minimum path length for the whole graph.}
  \item{Cp.rand}{mean of the clustering coefficient for the whole graph.}
  \item{in.degree.all}{vector of the degree of each node of the graph.}
  \item{Lp.rand.all}{vector of the minimum path length of each node of the graph.}
  \item{Cp.rand.all}{vector of the clustering coefficient of each node of the graph.}
}
\references{
S. H. Strogatz (2001)
Exploring complex networks. \emph{Nature}, Vol. 410, pages 268-276.

 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 

}
\author{S. Achard}
\seealso{\code{\link{rand.sw}},\code{\link{equadist.rand.sw}}  }
\examples{
mat<-sim.reg(8,20)
result<-reg.sw(8,20,dist.mat=matrix(1,8,8))
}
\keyword{graphs}

