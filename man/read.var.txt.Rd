\name{read.var.txt}
\alias{read.var.txt}
\alias{save.var.txt}
\title{Exportation and importation internal objects.}
\description{
Reads text files and imports them in object of class \code{"Wave Variance"}.
Exports object of class \code{"Wave Variance"}.
}
\usage{
read.var.txt()
save.var.txt(wave.var.list)
}
\arguments{
  \item{wave.var.list}{object of class \code{"Wave Variance"} containing the vectors of variance to be exported.}
}
\details{
These two functions cannot be used separately.
The names of the files used in \code{read.var.txt} are given by the \code{save.var.txt} functions
}
\value{
For \code{read.var.txt} function : Object of class \code{"Wave Variance"}, basically, a list with the following
  components
 \item{d?}{vectors of variance for each scale of the wavelet decomposition.}
 \item{lowerd?}{vector containing the lower bound of the variance for each scale of the wavelet decomposition.}
 \item{upperd?}{vector containing the upper bound of the variance for each scale of the wavelet decomposition.}
}
\references{}
\author{S. Achard}
\note{change between version 1 and 2, now the length of the time series and the number of the time series are saved with 
the values of the correlation.}
\seealso{\code{\link{read.cor.txt}}, \code{\link{save.cor.txt}}}
\examples{
data(brain)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]


# Construction of the correlation matrices for each level of the wavelet decomposition
wave.var.list<-const.var.list(brain, method = "modwt" ,wf = "la8", n.levels = 4,
                               boundary = "periodic", export.data=FALSE)

#Export the data
save.var.txt(wave.var.list)

#Import the data
read.var.txt()

}
\keyword{utilities}

