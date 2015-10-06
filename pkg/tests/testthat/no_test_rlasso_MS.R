# DGP
library(hdm2)
n <- 100; p <- 10

X <- matrix(rnorm(n*p), ncol=p)
y <- X%*%rep(c(1,0), c(5,5)) + rnorm(n)



output1 <- rlasso(X,as.vector(y))
output2 <- rlasso(X,as.matrix(y))
output3 <- rlasso(as.list(X),y)
output4 <- rlasso(as.data.frame(X),y)


output5 <- rlasso(X,y)
output6 <- rlasso(X,y, penalty=list(method = "standard", lambda.start = NULL, c = 1.1, gamma = 0.1))
uoutput7 <- rlasso(X,y, penalty=list(method = "standard"))
output8 <- rlasso(X,y, penalty=list(method = "X-dependent", numSim=5000))
output9 <- rlasso(X,y, penalty=list(method = "X-dependent", c=1.2))
output10 <- rlasso(X,y, penalty=list(method = "X-independent", c=1.2))
output11 <- rlasso(X,y, penalty=list(method = "none", lambda.start=rep(1,p)))
output12 <- rlasso(X,y, penalty=list(method = "CV"))

# check methods

t1 <- predict(output2)
