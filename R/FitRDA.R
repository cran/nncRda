FitRDA <-
function(X, y, Xt, yt, alpha=seq(0,0.99,0.11), delta=seq(0,3,0.2), ...){
#
# Finds Misclassification error rate for trian and test data using RDA. 
#
# X        : Training sample, expression matrix
# Xt       : Test sample, expression matrix
# y        : Class in training data. 
# yt       : Class in test data. 
# alpha    : alpha vector in function 'rda'
# delta    : delta vector in function 'rda'
#
    stopifnot(is.vector(y), is.vector(yt), 
                ncol(X)==length(y), 
                ncol(Xt)==length(yt),
                nrow(Xt)==nrow(X))
    fit <- rda(X, y)
    sink("Junk.txt")
    tryCatch(
            fit.cv <- rda.cv(fit,X,y,alpha=alpha,delta=delta,...),
            finally = {
                     sink()
                     unlink("Junk.txt")
                     })
    a<-fit.cv$cv.err
    b<-fit.cv$ngene
    La<-length(alpha)
    Ld<-length(delta)
    IndDel<-rep(delta,each=La)
    IndAlp<-rep(alpha,Ld)
    Ia<-1:length(a)
    NInd<-Ia[a==min(a)] # finding the (alpha, delta) that produce(s) minimum cv.err
    Ngene<-b[NInd]  # the number of genes corresponding to selected (alpha, delta)
    OptInd<-NInd[Ngene==min(Ngene)] # Index corresponding to optimal selection.
    Alpha<-IndAlp[OptInd]
    Delta<-IndDel[OptInd] 
    nAlpha<-length(Alpha)
    ytp<- sapply(1:nAlpha,function(i) predict(fit,x=X,y=y,xnew=Xt,alpha=Alpha[i],delta=Delta[i]))   
# remark: if test sample only has one row
    if(length(yt)==1) 
        {MisclassTest<- ytp!=yt}  
    else
        {MisclassTest<-apply(ytp,2,function(x) sum(yt!=x)/length(yt))} 
    c(Train=min(a)/length(y), Test=mean(MisclassTest)) 
 }
