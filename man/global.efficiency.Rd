\name{efficiency}
\alias{global.efficiency}
\alias{local.efficiency}
\alias{global.cost}
\alias{cost.evaluator}
\title{Graph efficiency}
\description{
Computes various measures of efficiency of a graph using the definition given by (Latora, 2001 and 2003)
}
\usage{
global.efficiency(adj.mat, weight.mat)
local.efficiency(adj.mat, weight.mat)
global.cost(adj.mat, weight.mat)
cost.evaluator(x)
}
\arguments{
\item{adj.mat}{adjacency matrix of the graph}
\item{weight.mat}{weighted matrix associated to the graph}
\item{x}{real}
}
\details{
Formula for the global efficiency :
\deqn{E_{global} = \frac{1}{N(N-1)} \sum_{i \neq j \in G} \frac{1}{L_{i,j}}}{%
E_global = 1/N(N-1) sum_{i <= j in G} 1/L_{i,j}}
where \eqn{L_{i,j}} is the minimum path length between each pair of nodes, and G the graph.

Formula for the nodal efficiency for the node \eqn{i}:
\deqn{E_{nodal}(i) = \frac{1}{(N-1)} \sum_{j \in G} \frac{1}{L_{i,j}}}{%
E_nodal(i) = 1/(N-1) sum_{j in G} 1/L_{i,j}}

Formula for the local efficiency for node \eqn{i}:
\deqn{E_{local} = \frac{1}{N_{G{_i}}(N_{G{_i}}-1)} \sum_{j,k \in G_i} \frac{1}{L_{j,k}}}{%
E_local = 1/N_{G{_i}}(N_{G{_i}}-1) sum_{j,k in G_i} 1/L_{j,k}}
where \eqn{G_i} is the set of nodes which are nearest-neighbours of the \eqn{i}th node, i.e., the nodes in \eqn{G_i} are each directly connected by a single edge to the index node \eqn{i} and \eqn{N_{G{_i}}} is the number of nodes in the sub-graph \eqn{G_i}. 

The computation of the cost requires the definition of an internal function, called \code{cost.evaluator}. For the moment, the \code{cost.evaluator} is the identity. Refer to Latora (2001) for 
the exact definition and usage of this function.
}
\value{
\item{global.efficiency$nodal.eff}{vector containing the nodal efficiency for each node of the graph
(equal to the inverse of the harmonic mean of the path length, when two nodes are disconnected, 
the path length is taken to be infinity so the inverse is 0)}
\item{global.efficiency$eff}{real corresponding to the mean of the nodal efficiency for the 
whole graph}
\item{local.efficiency$loc.eff}{vector containing the local efficiency for each node of the graph
(see details for the exact fomula)}
\item{global.efficiency$eff}{real corresponding to the mean of the local efficiency for the 
whole graph}
\item{global.cost}{real corresponding to the mean of the cost for the 
whole graph}
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
\note{only in version 2 and higher
}
\seealso{ \code{\link{const.adj.list}}, \code{\link{small.world}}}
\examples{

data(brain)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]




n.regions<-dim(brain)[2]

#Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 6, 
                               boundary = "periodic", p.corr = 0.975)


sup.seq<-((1:10)/10) #sequence of the correlation threshold 
nmax<-length(sup.seq)
Eglob<-matrix(0,6,nmax)
Eloc<-matrix(0,6,nmax)
Cost<-matrix(0,6,nmax)

n.levels<-6

#For each value of the correlation thrashold
for(i in 1:nmax){
n.sup<-sup.seq[i]

#Construction of the adjacency matrices associated to each level of the wavelet decomposition
wave.adj.list<-const.adj.list(wave.cor.list, sup = n.sup)


#For each level of the wavelet decomposition
for(j in 1:n.levels){

Eglob.brain<-global.efficiency(wave.adj.list[[j]],
			weight.mat=matrix(1,n.regions,n.regions))
Eglob[j,i]<-Eglob.brain$eff

Eloc.brain<-local.efficiency(wave.adj.list[[j]],
				weight.mat=matrix(1,n.regions,n.regions))
Eloc[j,i]<-Eloc.brain$eff

Cost.brain<-global.cost(wave.adj.list[[j]],
				weight.mat=matrix(1,n.regions,n.regions))
Cost[j,i]<-Cost.brain

}}


plot(sup.seq,(1:nmax)/2,type='n',xlab='Correlation threshold, R',ylab='Global efficiency',
     cex.axis=2,cex.lab=2,xlim=c(0,1),ylim=c(0,1))


for(i in 1:n.levels){
lines(sup.seq,Eglob[i,],type='l',col=i,lwd=2)
}
legend(x="topright",legend=c("Level 1","Level 2","Level 3","Level 4",
				"Level 5","Level 6"),fill=TRUE,col=(1:n.levels),lwd=2)


}
\keyword{graphs}