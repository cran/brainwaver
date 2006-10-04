\name{compute.FDR}
\alias{compute.FDR}
\title{ False Discovery Rate computation }
\description{
 Computation of the p-value cut-off which controls the false discovery rate when the test statistics have positive regression dependency on each of the test statistics corresponding to the true null hypotheses.
}
\usage{
compute.FDR(pvalue.vec, q)
}
\arguments{
  \item{pvalue.vec}{a vector containing the p-value for each hypothesis test.}
  \item{q}{value of the desired False Discovery Rate, exactly the upper limit for the expectation of the proportion of false positives.}
}
\details{
 This code implements the FDR procedure described in Benjamini and Yekutieli (2001).
}
\value{
  a real giving the p-value cutt-off.
}
\references{
Benjamini Y. and Yekutieli D. (2001)
The control of the false discovery rate in multiple testing under dependency. \emph{The Annals of Statistics}, Vol. 29, No. 4, pages 1165-1188
}
\author{S. Achard}
\note{The GeneTS package have also an implementation of this function }
\seealso{}
\examples{

data(young)
brain<-as.matrix(brain)

# WARNING : To process only the first five regions
brain<-brain[,1:5]


# Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 4, 
                               boundary = "periodic", p.corr = 0.975)

# For scale 4
pvalue.cor<-p.value.compute(wave.cor.list[[4]],proc.length=dim(brain)[1], 
                             sup=0.44, num.levels=4)

# Computation of the p-value threshold using FDR procedure
pvalue.thresh<-compute.FDR(pvalue.cor,0.05)
}
\keyword{multivariate}
\keyword{htest}
