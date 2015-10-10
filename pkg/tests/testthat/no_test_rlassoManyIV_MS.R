# DGP
#library(hdm2)
library(hdm)
data <- DGP.HCHIV()

y <- data$y
d <- data$d
Z <- data$Z
X <- data$X
output1 <- rlassoIV(X,d,y,Z, post=TRUE, intercept=TRUE)
coef(output1)
print(output1)
summary(output1)
confint(output1)

output2a <- rlassoIV(X,d,y,Z, select.X=TRUE, select.Z=FALSE)
coef(output2a)
print(output2a)
summary(output2a)
confint(output2a)

output2aa <- rlassoIVselectX(X,d,y,Z, select.X=TRUE, select.Z=FALSE)
coef(output2aa)
print(output2aa)
summary(output2aa)
confint(output2aa)

output2b <- rlassoIV(X,d,y,Z, select.X=TRUE, select.Z=FALSE,  post=TRUE, intercept=TRUE)
coef(output2b)
print(output2b)
summary(output2b)
confint(output2b)

output3 <- rlassoIV(X,d,y,Z, select.X=FALSE, select.Z=TRUE)
coef(output3)
print(output3)
summary(output3)
confint(output3)

output4 <- rlassoIV(X,Z,y,d, select.X=FALSE, select.Z=FALSE)
coef(output4)
print(output4)
summary(output4)
confint(output4)

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

output5 <- rlassoEffect(X,y, index=c(1,3,45))