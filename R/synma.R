synma <-function(n=c(100, 100), nt=c(500,500), rho=c(0.9,0.9), B=c(20,20), m=c(100,100), fE=0.02, E=0.5, Si=1){
stopifnot(length(n)==2, length(nt)==2, length(rho)==2, length(B)==2, length(m)==2)
stopifnot(length(fE)==1, length(E)==1, length(Si)==1, fE<1, fE>0, Si>0)
#Simulates One Microarray Data with n[1]+n[2] subjects, G genes.
# Expression matrix
#output is list with components X, y, Xt, yt, where
#  X training, G-by-(n[1]+n[2])
#  y training, (n[1], n[2])
#  Xt test, G-(nt[1]+nt[2])
#  y test, nt[1]+nt[2])
#
#input parameters:
#n, nt - number of subjects, c(#"healthy", #"diseased")
#        in training and test samples
# rho - correlation parameter for "healthy" and "diseased"
# m - number of blocks for "healthy" and "diseased"
# B - block size for "healthy" and "diseased"
#So total number of genes is B[1]*m[1]. If B[2]*m[2] is less than G,
# m[2] is increased so that B[2]*m[2] is greater than or equal
# to G. If B[2]*m[2] is greater than G, the first G rows & columns
# are selected.
# fE - fraction of expressed genes
# E - expression level of expressed genes
#So the number of expressed genes is B*m*fE, 40 with current defaults
# Si - variance inflation factor for expressed genes
#
#Default parameter settings replicate Example 4.3
#
####################################################
#Utility functions
#
is.CovarianceMatrix <-function(S)
    all(t(S)==S) & all(svd(S)$d > 0)
# 
#SigGenes - covariance matrix of gene expressions
SigGenes <-
function(rho, B=100, m=100){
# Gene covariance matrix given in Guo et. al (2006, eqn 4.1)
#  with default settings as in Example 4.3
# B: block size
# m: number of blocks
#So the number of genes is m*B
    SS <- array(0, dim=c(B,B,m))
    for (h in 1:m) {
        if (h%%2 == 1) 
            r = rho
        else
            r = -rho
        SS[,,h] <- toeplitz(r^(0:(B-1)))
        }
    S<-diag(0, m*B)
    for (h in 1:m)
        for (i in 1:B)
            for (j in 1:B) 
                S[(h-1)*B+i, (h-1)*B+j] <- SS[i,j,h]
    S 
}
#####################################################
    G <- B[1]*m[1] #number of genes
    EG <- ceiling(fE*G) #number of expressed genes
    mu1 <- rep(0, G)
    mu2<- c(rep(E, EG), rep(0, (G-EG))) 
    SIG1 <- SigGenes(rho=rho[1], B=B[1], m=m[1])
    if (!is.CovarianceMatrix(SIG1))
        stop("SIG2 not valid covariance matrix!") 
    Smul <- c(rep(Si,EG), rep(1,G-EG))
    SIG2 <-((SigGenes(rho=rho[2], B=B[2], m=m[2]))[1:G[1],1:G[1]])*outer(Smul, Smul)
    if (!is.CovarianceMatrix(SIG2))
        stop("SIG2 not valid covariance matrix!")
    X1 <- mvrnorm(n[1], mu1, SIG1) 
    X2 <- mvrnorm(n[2], mu2, SIG2)
    X <- rbind(X1, X2)
    y <- rep(1:2, n)
    Xt1<-mvrnorm(nt[1], mu1, SIG1)
    Xt2<-mvrnorm(nt[2], mu2, SIG2)
    Xt<-rbind(Xt1, Xt2)
    yt <- rep(1:2, nt)
    list(X=t(X), y=y, Xt=t(Xt), yt=yt)
}
