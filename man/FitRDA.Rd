\name{FitRDA}
\alias{FitRDA}
\title{ Regularized Discriminant Analysis }
\description{
Fits RDA, computes mis-classification rates on training and test data.
}
\usage{
FitRDA(X, y, Xt, yt, alpha = seq(0, 0.99, 0.11), delta = seq(0, 3, 0.2), ...)
}

\arguments{
  \item{X}{training, expression matrix }
  \item{y}{response vector taking values 1, 2, ...}
  \item{Xt}{test, expression matrix }
  \item{yt}{ response vector for test data taking values 1, 2, ... }
  \item{alpha}{ tuning parameter }
  \item{delta}{ turning parameter }
  \item{...}{ optional arguments, passed to rda }
}

\value{
Error rates for training and test data
}


\references{ 
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\examples{
#Tiny test example - 10 seconds
M<-synma(n=c(10,10), nt=c(10,10), B=c(5,5), m=c(5,5))
FitRDA(M$X, M$y, M$Xt, M$yt)
}
\keyword{ classif }
