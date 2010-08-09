\name{read.cor.txt}
\alias{read.cor.txt}
\alias{save.cor.txt}
\title{Exportation and importation of internal objects.}
\description{
Reads text files and imports them in 
object of class \code{"Wave Correlation"}.
Exports object of class \code{"Wave Correlation"}.
}
\usage{
read.cor.txt()
save.cor.txt(wave.cor.list)
}
\arguments{
  \item{wave.cor.list}{object of class \code{"Wave Correlation"} containing the correlation matrices to be exported.}
}
\details{
These two functions cannot be used separately.
The names of the files used in \code{read.cor.txt} are given by the \code{save.cor.txt} functions
}
\value{
For \code{read.cor.txt} function : Object of class \code{"Wave Correlation"}, basically, a list with the following
  components
 \item{d?}{Correlation matrix for each scale of the wavelet decomposition.}
 \item{lowerd?}{matrix containing the lower bound of the correlation for each scale of the wavelet decomposition.}
 \item{upperd?}{matrix containing the upper bound of the correlation for each scale of the wavelet decomposition.}
}
\author{S. Achard}
\note{change between version 1 and 2, now the length of the time series and the number of time series are saved with 
the values of the correlation.}
\seealso{\code{\link{read.var.txt}}, \code{\link{save.var.txt}}}
\examples{
data(brain)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]



# Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 4, 
                              boundary = "periodic", p.corr = 0.975,export.data=FALSE)

#Export the data
save.cor.txt(wave.cor.list)

#Import the data
read.cor.txt()

}
\keyword{utilities}

