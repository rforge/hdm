set.seed(1)
n = 100 #sample size
p = 100 # number of variables
s = 3 # nubmer of variables with non-zero coefficients
X = Xnames = matrix(rnorm(n*p), ncol=p)
colnames(Xnames) <- paste("V", 1:p, sep="")
beta = c(rep(3,s), rep(0,p-s))
Y = X%*%beta + rnorm(n)
Xnew = Xnewnames = matrix(rnorm(n*p), ncol=p)  # new X
colnames(Xnewnames) <- paste("V", 1:p, sep="")
Ynew =  Xnew%*%beta + rnorm(n)  #new Y
dat <- cbind(Y,Xnames)
colnames(dat)[1] <- "Y"
dat <- as.data.frame(dat)

Xd <- X[,-1]
d <- X[,1]
reg1 <- rlassoEffect(Xd,Y,d)
summary(reg1)
confint(reg1)
reg2 <- rlassoEffect(X[,-1],Y,X[,1])
library(ggplot2)
plot(reg1)
##########################################################################

from1 <-  paste(paste("V", 1:p, sep=""), collapse="+")
from2 <- paste("Y ~" , from1)
reg3 <- rlassoEffects(X,Y, index=c(1,2,5))
plot(reg3)
reg4 <- rlassoEffects(as.formula(from2), data = dat, I = ~ V1 + V4 + V30)
