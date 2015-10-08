# DGP
library(hdm)
n <- 100; p <- 10

X <- matrix(rnorm(n*p), ncol=p)
y <- X%*%rep(c(1,0), c(5,5)) + rnorm(n)



output1 <- rlassoLM(X,y, index=1:10)
print(output1)
summary(output1)
confint(output1, joint=FALSE)
confint(output1, joint=TRUE)
plot(output1)


 output2 <- rlasso(X,as.matrix(y))
output3 <- rlasso(as.list(X),y)
output4 <- rlasso(as.data.frame(X),y)

# check methods

