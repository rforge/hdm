# DGP
library(hdm2)
data <- DGP.IV()

y <- data$y
d <- data$X
Z <- data$Z
X <- NULL
output1 <- rlassoIV(X,d,y,Z, post=TRUE, intercept=TRUE)
coef(output1)
print(output1)
summary(output1)
confint(output1)

output2 <- rlasso(X,as.matrix(y))
output3 <- rlasso(as.list(X),y)
output4 <- rlasso(as.data.frame(X),y)


output5 <- rlasso(X,y)
output6 <- rlasso(X,y, penalty=list(method = "standard", lambda.start = NULL, c = 1.1, gamma = 0.1))
output7 <- rlasso(X,y, penalty=list(method = "standard"))
output8 <- rlasso(X,y, penalty=list(method = "X-dependent", numSim=5000))
output9 <- rlasso(X,y, penalty=list(method = "X-dependent", c=1.2))
output10 <- rlasso(X,y, penalty=list(method = "X-independent", c=1.2))
output11 <- rlasso(X,y, penalty=list(method = "none", lambda.start=rep(1,p)))
output12 <- rlasso(X,y, penalty=list(method = "CV"))

# check methods

