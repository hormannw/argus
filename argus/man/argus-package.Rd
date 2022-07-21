\encoding{UTF-8}
\name{argus-package}
\alias{argus-package}
\alias{argus}
\docType{package}

\title{
  Generator and density for the Argus distribution
}

\description{
  This package provides two random variate generation methods, the density, the CDF and the quantile function
of the argus distribution.
}

\details{
  \tabular{ll}{
    Package: \tab argus          \cr
    License: \tab GPL 2 or later\cr
  }

  Package \pkg{argus} provides the routines:
  \describe{
    \item{\code{\link{rargus}}}{generates argus distributed random variates.}
    \item{\code{\link{dargus}}}{computes the density of the argus distribution.}
    \item{\code{\link{pargus}}}{computes the CDF of the argus distribution.}
    \item{\code{\link{qargus}}}{computes the quantile function (ie. the inverse CDF )of the argus distribution.}
  }

}

\author{
  Wolfgang \enc{Hörmann}{Hoermann} & Christoph Baumgarten
}

\references{%
 C. Baumgarten, W. Hörmann: Generating ARGUS random variates
}

\examples{
## Evaluate the pdf (density)
dargus(c(0.1,0.5,0.9), chi=0.3)
# chi must have the same length as x !
dargus(c(0.1,0.5,0.9), chi=c(0.3,1.3,2.3)) 
## Evaluate the CDF
pargus(c(0.1,0.5,0.9), chi=0.3)
pargus(c(0.1,0.5,0.9), chi=c(0.3,1.3,2.3)) 
## Evaluate the quantile function
pargus(c(0.1,0.5,0.9), chi=4.5)
## Draw a random sample
rargus(n=10, chi=0.3)
## compare histogram and density
system.time(y<-rargus(1.e6,chi=2.5))
hist(y,breaks=100,freq=FALSE)
lines(x<-seq(0,1,1.e-3),dargus(x,2.5),col=2)
}

\keyword{package}
\keyword{distribution}
\keyword{datagen}
