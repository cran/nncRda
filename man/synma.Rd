\name{synma}
\alias{synma}
\title{
Generate Synthetic Microarray Data
}
\description{
Simulates data from a microarray
}
\usage{
synma(n = c(100, 100), nt = c(500, 500), rho = c(0.9, 0.9), B = c(20, 20), m = c(100, 100), fE = 0.02, E = 0.5, Si = 1)
}

\arguments{
  \item{n}{
vector length 2, c(number healthy, number diseases), training sample
}
  \item{nt}{
vector length 2, c(number healthy, number diseases), test sample
}
  \item{rho}{
vector length 2, c(correlation healthy, correlation diseases)
}
  \item{B}{
vector length 2, c(block size healthy, block size disease)
}
  \item{m}{
vector length 2, c(# blocks healthy, # blocks disease)
}
  \item{fE}{
fraction of genes expressed in 'diseased' patients
}
  \item{E}{
Expression level of diseased genes
}
  \item{Si}{
Variance inflation factor applied to the 'diseased' group
}

}
\details{
This is a generalization of the simulation model for microarray data
suggested by Guo et al. (2007).
We have extended the simulation to allow for different covariance
structures across groups.
This occurs with actual microarray data as in
the lymphoma data discussed by
Dudoit et al. (2002, Section 5.1.3 and Table 1).
}

\value{
A list with four elements:

\item{X}{expression matrix, training, (B*m)-by-(n[1]+n[2]), where B=B[1] and m=m[1]}
\item{y}{response, training, n[1] healthy followed by n[2] diseased}
\item{Xt}{expression matrix, test, (B*m)-by-(nt[1]+nt[2]), where B=B[1] and m=m[1]}
\item{yt}{response, test, nt[1] healthy followed by nt[2] diseased}
}

\references{ 
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 

Dudoit, S., Fridlyand, J. and Speed, T. (2002). Comparison of discrimination
methods for the classification of tumors using gene expression data. J. Amer. Stat. Assoc.,
97, 7787.


Islam, M. S. (2008). 
Periodicity, Change Detection and Prediction in Microarrays. Ph.D. thesis, 
The University of Western Ontario, Canada.
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\seealso{\code{\link[nnc]{nnc}}}
\examples{
#Example 1. A micro-sized microarray!
ans<-synma(n=c(3,3),nt=c(2,2),B=c(3,3),m=c(3,3),fE=0.5)
ans

#Example 2. Synethetic data 
\dontrun{
#Source: TestA.R
#about 140 sec
#
#40*50 = 2000 genes
#training sample size 200 (100 in each class)
#test sample 1000 (500 in each class)
#fE=0.05 mean 5% of genes are expressed at 0.5 rest at 0.0
#Almost the same as in thesis, Table 3.1
#kNN autocovariate, optimal k vs k=1 with k=0 too 
set.seed(321777)
StartTime <- proc.time()[3]
M<-synma(n=c(100,100), nt=c(500,500), B=c(40,80), m=c(20,10), fE=0.05)
ANS1<-nncRda(M$X, M$y, M$Xt, M$yt)
ANS2<-nncRda(M$X, M$y, M$Xt, M$yt, k=1)
ANS<-c(ANS1,ANS2[2])
names(ANS)<-c(names(ANS1), names(ANS2)[2])
EndTime <- proc.time()[3]
EndTime-StartTime
ANS
}
}
\keyword{ classif }

