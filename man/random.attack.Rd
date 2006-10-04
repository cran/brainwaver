\name{graph.attack}
\alias{random.attack}
\alias{node.attack}
\alias{targeted.attack}
\title{Attack of graphs}
\description{
Computes the evolution of the size of the largest connex component and its mean minimum path length under a random or targeted attack of a graph.
}
\usage{
random.attack(adj.mat, max.remove, nsim)
node.attack(adj.mat, node.sup)
targeted.attack(adj.mat, max.remove)

}
\arguments{
  \item{adj.mat}{adjacency matrix of the graph.}
  \item{max.remove}{maximum number of nodes to be removed}
  \item{nsim}{number of simulation of random attack to be processed}
\item{node.sup}{number of the node to be removed}
}
\details{
An attack consists in removing a node and all its connections from a given graph. \code{random.attack} removes \code{max.remove} nodes by selecting one regional node at random and removing it (and all its connections) from the graph. \code{targeted.attack} removes \code{max.remove} nodes : the first node to be eliminated was the hub with the largest degree, and nodes were subsequently eliminated in rank order of decreasing degree.
}
\value{
\item{size.large.connex}{vector containing the evolution of the size of the largest connexe component during an attack.}
\item{charac.path.length}{vector containing the evolution of the mean minimum path length of the largest connexe component during an attack.}
}
\references{
Albert R., Barabasi A. L. (2002) 
Statistical mechanics of complex networks. 
\emph{Rev Mod Physics} 74 pages 47-97. 

 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
}
\author{S. Achard}
\note{}
\seealso{}
\examples{

set2<-array(c(5,6.5,7,6.5,5,3.5,3,3.5,1,1.5,3,4.5,5,4.5,3,1.5),dim=c(8,2))
names<-c(1:8)

mat<-sim.rand(8,20)

#simulated attack on node 1
plot(set2[,1], set2[,2], type = "n",xlab="", ylab="",cex.lab=1.5)
text(set2[2:8,1], set2[2:8,2], names[2:8], cex = 1.5)
text(set2[1,1], set2[1,2], 1 , cex = 1.5,col="gray")

for(k in 2:8){
                for(q in 1:(k-1)){

                if(mat[k,q]==1)
                {
                if(q==1) visu <- "gray" else visu <- "red"
                lines(c(set2[k,1], set2[q,1]), c(set2[k,2], set2[q,2]), col = visu)
        }

}
}

ra<-random.attack(mat,8,10)
na<-node.attack(mat, 1)
ta<-targeted.attack(mat,8)


}
\keyword{graphs}

