\name{nncRda}
\alias{nncRda}
\title{Test error rate for rda with nnc}
\description{This finds error rate for rda and improved rda with nnc. If kSelect=TRUE, then the optimal nearest neighbour covariate is used; otherwise, only nnc corresponding to first nearest neighbour is used.}
\usage{
nncRda(X, y, Xt, yt, k = "default", alpha = seq(0, 0.99, 0.11), delta = seq(0, 3, 0.2), ...)
}

\arguments{
  \item{X}{training, expression matrix}
  \item{y}{training output}
  \item{Xt}{test expression matrix}
  \item{yt}{test output}
  \item{k}{"default" use MLE k. Otherwise set to specified value, k>0.}
  \item{alpha}{ alpha in rda function}
  \item{delta}{ delta in rda function}
  \item{\dots}{ optional arguments passed to rda}
}
\details{
cross validation
}

\value{A vector representing error rate for rda and improved rda with nnc.}
\references{
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 

Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays. Ph.D. thesis, The University of Western Ontario, Canada.
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\seealso{\code{\link[nnc]{nnc}}}

\examples{
#Example 1.
#Tiny test example - 10 seconds
M<-synma(n=c(10,10), nt=c(10,10), B=c(5,5), m=c(5,5))
nncRda(M$X, M$y, M$Xt, M$yt)

#Example 2. Takes about 120 sec
\dontrun{
#40*20 = 800 genes
#training sample size 200 (100 in each class)
#test sample 1000 (500 in each class)
#fE=0.05 mean 5% of genes are expressed at 0.5 rest at 0.0
#kNN autocovariate, optimal k vs k=1 with k=0 too 
set.seed(321777)
StartTime <- proc.time()[3]
M<-synma(n=c(100,100), nt=c(500,500), B=c(40,80), m=c(20,10), fE=0.05)
ANS1<-nncRda(M$X, M$y, M$Xt, M$yt)
ANS2<-nncRda(M$X, M$y, M$Xt, M$yt, k=1)
ANS<-c(ANS1,ANS2[2])
names(ANS)<-c(names(ANS1), names(ANS2)[2])
EndTime <- proc.time()[3]
TotalTime<-EndTime-StartTime
TotalTime
ANS
sqrt(ANS*(1-ANS)/1000)
##elapsed 
## 123.21 
##> ANS
##    RDA RDA,k=5 RDA,k=1 
##  0.363   0.071   0.092 
##> sqrt(ANS*(1-ANS)/1000)
##        RDA     RDA,k=5     RDA,k=1 
##0.015206282 0.008121515 0.009139803 
##Conclusion: the difference in misclassification rates is 
##   quite large and highly signficant
}

}

\keyword{ classif }


