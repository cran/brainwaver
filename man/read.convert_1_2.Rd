\name{read.convert_1_2}
\alias{read.convert_1_2}
\title{Conversion of already edited objects from version 1 to 2}
\description{
Reads text files (\code{"wave_cor_mat_level_[1:n.levels].txt"}) from version 1 and 
add the length of the 
time series in the object of class \code{"Wave Correlation"}.
}
\usage{
read.convert_1_2(proc.length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{proc.length}{ specifies the length of the original processes using to construct the \code{cor.mat}}~~Describe \code{proc.length}
}
\details{
}
\value{
 Object of class \code{"Wave Correlation"}, basically, a list with the following
  components
 \item{d?}{Correlation matrix for each scale of the wavelet decomposition.}
 \item{lowerd?}{matrix containing the lower bound of the correlation for each scale of the wavelet decomposition.}
 \item{upperd?}{matrix containing the upper bound of the correlation for each scale of the wavelet decomposition.}
}
\references{}
\author{S. Achard}
\note{only in version 1.2}
\seealso{\code{\link{read.cor.txt}}, \code{\link{save.cor.txt}}}
\keyword{utilities}

