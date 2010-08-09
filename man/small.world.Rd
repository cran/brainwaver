\name{small.world}
\alias{small.world}
\title{Computation of some quantitive parameters given in the framework of the small-world graph theory.}
\description{
For a given graph, computes the size of the largest connex component, degree, minimum path length and clustering coefficient.
}
\usage{
small.world(wave.adj.mat, dat = "reduced", distance = "norm", 
            coord = 0, export.data = FALSE)
}
\arguments{
  \item{wave.adj.mat}{ adjacency matrix of the graph to be analysed.}
  \item{dat}{ character, if \code{dat = "all"}, the degree, minimum path length and clustering coefficient are computed for each node of the graph. If \code{dat = "reduced"}, only the mean on all the nodes is computed.}
  \item{distance}{matrix with a distance associated to each pair of nodes of the graph to take into account in the computation of the small-world parameters. By default, the matrix of the distance has a value of 1 in each entry.}
  \item{coord}{optional vector containing the coordinate of the nodes of the graph. }
 \item{export.data}{logical. If TRUE the correlation matrices with the upper and lower bound are exported to text file.}
}
\value{
  \item{in.degree.mean}{ real corresponding to the mean of the degree of all the nodes of the graph.}
 \item{in.degree}{vector containing the degree of each node of the graph.}
\item{Lp.mean}{ real corresponding to the mean of minimum path length for each node belonging to the largest connex component of the graph only.}
 \item{Lp}{vector containing the mean minimum path length for each node of the graph. (the computaion of Lp requires at least one neighbourgh per node, if it is not the case the value of Lp is equal to -1)}
\item{Cp.mean}{real corresponding to the mean of clustering coefficient for each node belonging to the largest connex component of the graph only.}
 \item{Cp}{vector containing the clustering coefficient for each node of the graph. (the computaion of Cp requires at least two neighbourghs per node, if it is not the case the value of Cp is equal to -1)}  
\item{size.large.connex}{real corresponding to the size of the largest connex component.}
}
\references{
D. J. Watts and S. H. Strogatz (1998)
Collective dynamics of ``small-world'' networks. \emph{Nature}, Vol. 393, pages 440-442.

S. H. Strogatz (2001)
Exploring complex networks. \emph{Nature}, Vol. 410, pages 268-276.

 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 



}
\author{S. Achard}
\examples{
# fig 3 of Achard (2006)


data(brain)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]



#Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 6, 
                               boundary = "periodic", p.corr = 0.975)

sup.seq<-((1:10)/10) #sequence of the correlation threshold 
nmax<-length(sup.seq)
in.degree.mean<-matrix(0,6,nmax)
n.levels<-6

#For each value of the correlation thrashold
for(i in 1:nmax){
n.sup<-sup.seq[i]

#Construction of the adjacency matrices associated to each level of the wavelet decomposition
wave.adj.list<-const.adj.list(wave.cor.list, sup = n.sup, proc.length=dim(brain)[1])


#For each level of the wavelet decomposition
for(j in 1:n.levels){

param.sw.brain<-small.world(wave.adj.list[[j]],dat="reduced")
in.degree.mean[j,i]<-param.sw.brain$in.degree.mean

}}

#Plots of the average in-degree in terms of the scale

n.regions<-dim(brain)[2]

plot(sup.seq,(1:nmax)/2,type='n',xlab='Correlation threshold, R',ylab='Mean degree, k',
     cex.axis=2,cex.lab=2,xlim=c(0,1),ylim=c(0,90))


for(i in 1:n.levels){
lines(sup.seq,in.degree.mean[i,],type='l',col=i,lwd=2)
}
lines(sup.seq,rep(log(n.regions),nmax))
legend(x="topright",legend=c("Level 1","Level 2","Level 3","Level 4",
				"Level 5","Level 6"),col=(1:n.levels),lwd=2)





}
\keyword{graphs}



