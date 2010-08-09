\name{rand.eff}
\alias{rand.eff}
\alias{reg.eff}
\title{Efficiency for simulated graphs}
\description{
Computes the local, global efficiency and cost for simulated 
random and regular graphs.
}
\usage{
rand.eff(nsim, n.nodes.rand, n.edges.rand, dist.mat, dat = "reduced")
reg.eff(n.nodes.rand, n.edges.rand, dist.mat)
}
\arguments{
\item{nsim}{number of simulated graphs to use for the computation of the small-world parameters.}
  \item{dat}{character string specifying if all the small-world parameters have to be returned. If \code{"reduced"}, only the mean of the parameters for the whole graph is returned.}
  \item{n.nodes.rand}{number of nodes of the simulated graphs}
  \item{n.edges.rand}{number of edges of the simulated graphs}
\item{dist.mat}{matrix with a distance associated to each pair of nodes of the graph to take into account in the computation of the efficiency values.}
}
\value{
\item{eff}{global efficiency for the whole graph}
\item{loc}{local efficiency for th whole graph}
\item{cost}{cost for th whole graph}
}
\references{ 
V. Latora, M. Marchiori (2001) 
Efficient Behavior of Small-World Networks. \emph{Phys. Rev. Lett.}, Vol. 87, N. 19, pages 1-4.

V. Latora, and M. Marchiori (2003) 
Economic Small-World Behavior in Weighted Networks. \emph{Europ. Phys. Journ. B}, Vol. 32, pages 249-263. 

S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
}
\author{S. Achard}
\note{only in version 1.2 and after}
\seealso{ \code{\link{const.adj.list}}, \code{\link{small.world}}}
\examples{

result<-rand.eff(10,8,20,dist.mat=matrix(1,8,8))
result<-reg.eff(8,20,dist.mat=matrix(1,8,8))

}
\keyword{graphs}