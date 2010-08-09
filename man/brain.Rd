\name{brain}
\alias{brain}
\docType{data}
\title{Time series obtained by an fMRI experiement on the brain}
\description{
  Time series for each region of interest in the brain. These series are obtained by SPM preprocessing.
}
\usage{data(brain)}
\format{
  A data frame with 2048 observations on the following 90 variables.
}
\source{
contact S. Achard (sa428@cam.ac.uk)
}
\references{
R. Salvador, J. Suckling, C. Shwarzbauer, E.T. Bullomre (2005)
Undirected graphs of frequency-dependent functional connectivity in the whole brain networks. \emph{Philos Trans R Soc Lond B Biol Sci} Vol. 360, pages 937-946.

 S. Achard, R. Salvador, B. Whitcher, J. Suckling, Ed Bullmore (2006)
A Resilient, Low-Frequency, Small-World Human Brain Functional Network with Highly Connected Association Cortical Hubs. \emph{Journal of Neuroscience}, Vol. 26, N. 1, pages 63-72. 


 }
\examples{
data(brain)
## maybe str(brain) ; plot(brain) ...
}
\keyword{datasets}
