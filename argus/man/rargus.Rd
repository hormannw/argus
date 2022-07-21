\encoding{UTF-8}

\name{rargus}
\alias{rargus}
\alias{dargus}
\alias{pargus}
\alias{qargus}

\title{ 
  The Argus distribution
}

\description{
Density, distribution function, quantile function and random generation for the Argus distribution with parameter \code{chi}.
\code{rargus()} uses very efficient random variate generation methods. Inversion or the Ratio of Uniforms method can be selected. 
}

\usage{
rargus(n=length(chi) , chi, method = c("inversion","RoU") )
dargus(x, chi, log = FALSE)
pargus(x, chi, lower = TRUE, log.p = FALSE)
qargus(p, chi, lower = TRUE, log.p = FALSE)}

\arguments{
  \item{n}{Number of observations}
  \item{chi}{parameter vector \code{chi>0}}
  \item{method}{random variate generation method to be used}
  \item{x}{x-value for the density or CDF}
  \item{log}{If \code{TRUE} the logarithm of the density is returned.}
  \item{lower}{If \code{FALSE} 1-CDF or its inverse is returned}
  \item{log.p}{\code{log(p)} the logarithm of the probability is used or returned}
  \item{p}{probability of the quantile}
}

\details{
  The generators used for \code{rargus()} are very efficient also for the varying parameter case,
  i.e. when the vector \code{chi} has length \code{n}.
  
  \code{method = "inversion"}: transforms the uniform variate into an argus variate
   using a close to exact approximation of the inverse CDF. It is also well suited for the 
   varying parameter case as the theorem that an argus variate can be represented as a transformed 
   truncated Gamma(1.5) variate allows to express the inverse CDF of the argus distribution as a 
   simple transform of the inverse CDF of the Gamma(1.5) distribution.
   The use of the \code{pinv.new()} function of the \code{Runuran} package makes that evaluation very fast. 
   When loading the argus-package \code{pinv.new()} calculates and stores the required tables.
   
   \code{method = "RoU"}: uses the Ratio of Uniforms method which requires 2 or more uniform variates
   to generate one argus variate. It is also a bit slower than the inversion method and is added here
   mainly for the case that the \code{Runuran}-package is not available.
}

\value{
\code{rargus} creates a random sample of size \code{n}. 
 
\code{dargus} gives the density, \code{dargus} give the CDF and \code{qargus} gives the quantile function.
The length of the result for these 3 functions is the maximum of the lengths of the numerical arguments that
are recycled to the length of the result. Only the first elements of the logical arguments are used.
}

\author{
  Wolfgang \enc{Hörmann}{Hoermann} and
  Christoph Baumgarten
}

\references{%
 C. Baumgarten, W. Hörmann: Generating ARGUS random variates
}

%% \seealso{}

\examples{
## Evaluate the pdf (density)
dargus(c(0.1,0.5,0.9), chi=c(0.3,1.3,2.3)) 
## Evaluate the CDF
pargus(c(0.1,0.5,0.9), chi=c(0.3,1.3,2.3)) 
## Draw a random sample
rargus(n=5, chi=0.3)
## Generate a large sample with a fixed parameter
system.time(y<-rargus(n=1.e6,chi=2.5))
## compare histogram and density
hist(y,breaks=100,freq=FALSE)
lines(x<-seq(0,1,1.e-3),dargus(x,chi=2.5),col=2,lwd=2)
## transform the sample with the CDF and check if it is then uniform
uy <- pargus(y,chi=2.5)
hist(uy,breaks=100)
## to demonstrate the speed also of the varying parameter case:
chiv<-3*runif(2.e6)
system.time(y<-rargus(n=2.e6,chi=chiv))
hist(pargus(y,chiv),breaks=100)

#\dontrun{
#}

}

\keyword{package}
\keyword{distribution}
\keyword{datagen}
