nncRda <-
function(X, y, Xt, yt, k="default", alpha=seq(0,0.99,0.11), delta=seq(0,3,0.2), ...){
#
#Misclassification rate for train and test data using RDA
#
# X        : Training sample, expression matrix
# Xt       : Test sample, expression matrix
# y        : Class in training data: y should be 1, 2...
# yt       : Class in test data. 
# k        : "default",  use pseudolikelihood, k>0 set k in kNN
# alpha    : alpha vector in function 'rda'
# delta    : delta vector in function 'rda'
#
	stopifnot(is.vector(y), is.vector(yt), 
				ncol(X)==length(y), 
				ncol(Xt)==length(yt),
	 			nrow(Xt)==nrow(X),
				is.character(k)|is.numeric(k))
	is.wholenumber <-
		function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
	if (is.numeric(k))
		stopifnot(is.wholenumber(k), k>=1)
	if(is.character(k)) 
		K <- khat(t(X), y, plot=FALSE) 
	else
		K <- k 
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
# find ALL the (alpha, delta) that produce(s) minimum cv.err
	Ia <- 1:length(a)
	NInd<-Ia[a==min(a)] 
	Ngene<-b[NInd]  
	OptInd<-NInd[Ngene==min(Ngene)]
	Alpha<-IndAlp[OptInd]
	Delta<-IndDel[OptInd] 
	nAlpha<-length(Alpha)
	ytp<- sapply(1:nAlpha, function(i) 
				predict(fit,x=X,y=y,xnew=Xt,alpha=Alpha[i],delta=Delta[i])
				)   
	if(length(yt)==1) 
		{MisclassTest<- ytp!=yt}  
	else
		{MisclassTest<-apply(ytp,2,function(x) sum(yt!=x)/length(yt))} # we do so in case test sample has one row
	TestRDA <- mean(MisclassTest) 
	XZ <- rbind(X, nnc(t(X),y, K))
	XZt <- rbind(Xt, nncTest(t(X), y, t(Xt), K))
	Testnnc <- FitRDA(XZ, y, XZt, yt)[2]
	res <- c(TestRDA, Testnnc)
	names(res) <- c("RDA", paste("RDA,k=", K, sep=""))
	res
	}
	
