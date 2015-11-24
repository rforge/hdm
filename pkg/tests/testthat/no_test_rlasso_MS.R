# DGP
library(hdm2)
n <- 100; p <- 10

X <- matrix(rnorm(n*p), ncol=p)
y <- X%*%rep(c(1,0), c(5,5)) + rnorm(n)



output1 <- rlasso(X,as.vector(y))
output2 <- rlasso(X,as.matrix(y))
output3 <- rlasso(as.list(X),y) # does not work
output4 <- rlasso(as.data.frame(X),y)


output5 <- rlasso(X,y)
output6a <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=FALSE, X.design="independent"))
output6b <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=FALSE, X.design="independent", lambda.start=rep(10,p)))
output7a <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=TRUE, X.design="independent"))
output7b <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=TRUE, X.design="independent", lambda.start=rep(10,p)))
output8a <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=FALSE, X.design="dependent"))
output8b <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=FALSE, X.design="dependent", lambda.start=rep(10,p)))
output9a <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=TRUE, X.design="dependent"))
output9b <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic=TRUE, X.design="dependent", lambda.start=rep(10,p)))
output10 <- rlasso(X,y, post=TRUE, intercept=TRUE, normalize=TRUE, penalty=list(homoscedastic="none", lambda.start=rep(10,p)))
# check methods

t1 <- predict(output2)
