\name{correlations.to.adjacencies}
\alias{correlations.to.adjacencies}
\alias{ideal.wavelet.levels}
\alias{distance}
\title{Produce adjencency matrices for a given number of edges}
\description{
 Given a correlations thingy as produced by \code{const.cor.list},
produce a list of adjacency matrices fiddled to have a preferred number of edges
Actually this is not quite possible, but come as close as \code{choose.thresh.nbedges}
will allow us. A functional parameter allows us to say things like produce the graphs with n log n edges where n is the number of nodes
}
\usage{
correlations.to.adjacencies(correlations, edge.func)
ideal.wavelet.levels(brain)
distance(x,y,z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{correlations}{a list of correlation matrices produced by \code{const.cor.list}}
  \item{edge.func}{a function to mention the way to choose the number of edges given the number of nodes in the graph. In the companion scripts files, the small-limit is used and by default \code{edge.func=(function(x){x*log(x)})}}
\item{brain}{ matrix containing the data time series. Each column of the matrix represents one time series.}
\item{x}{x coordinate}
\item{y}{y coordinate}
\item{z}{z coordinate}
}
\details{
Functions produced to manipulate better nice outputs of the package
}
\value{
  \item{correlations.to.adjacencies}{Description of 'comp1'}
   \item{ideal.wavelets.levels}{number indicating up to each wavelet scale it is possible to go given the length of the time series}
  \item{disctance}{the euclidean distance in 3D} 
}
\references{ S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 
 }
\author{John Aspden, external collaborator of the brainwaver package}
\note{
}
\seealso{\code{\link{const.cor.list}}}
\examples{

data(brain)
brain<-as.matrix(brain)
# WARNING : To process only the first five regions
brain<-brain[,1:5]

n.levels<-4
wave.cor.list<-const.cor.list(brain,n.levels=n.levels)
adj.mat<-correlations.to.adjacencies(wave.cor.list,edge.func=(function(x){x*log(x)}))


}
\keyword{multivariate}
\keyword{ts}
