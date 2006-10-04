\name{equadist.rand.sw}
\alias{equadist.rand.sw}
\title{Small-world parameters for a simulated graph given its degree distribution}
\description{
Computes the degree, the minimum path length and the clustering coefficient for a simulated graph with a known degree distribution.
}
\usage{
equadist.rand.sw(nsim, dat = "reduced", dist.mat, degree.dist)
}
\arguments{
  \item{nsim}{number of simulated graphs to use for the computation of the small-world parameters.}
  \item{dat}{character string specifying if all the small-world parameters have to be returned. If \code{"reduced"}, only the mean of the parameters for the whole graph is returned.}
  \item{dist.mat}{matrix with a distance associated to each pair of nodes of the graph to take into account in the computation of the small-world parameters.}
  \item{degree.dist}{vector describing the degree distribution verified by the nodes of the simulated graph.}
}
\details{
Because of the problem for the simulation of a random graph with exactly the same degree distribution, each simulation applied in this function is subject to the result to get a graph with the same number of edges as expected. So, the number of wanted simulations can be very different from the number of simulations taken into account in the computation of the small-world parameters.
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
\note{}
\seealso{\code{\link{sim.equadist}}, \code{\link{rand.sw}}, \code{\link{reg.sw}} }
\examples{

#For a scale-free graph

x<-1:50
probx<-x^(-1.4)
n.nodes<-50
n.edges<-250

stop<-0

while(stop==0){
write.table(stop)
r<-sample(x,n.nodes,prob=probx,replace=TRUE)
if(sum(r)==n.edges) stop<-1
  }

sf.degree<-r

mat<-sim.equadist(sf.degree)

result<-equadist.rand.sw(10, dat = "reduced", dist.mat=matrix(1,50,50), 
	                  degree.dist=sf.degree)

}
\keyword{graphs}


