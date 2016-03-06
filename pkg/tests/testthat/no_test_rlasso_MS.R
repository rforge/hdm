set.seed(1)
n = 100 #sample size
p = 100 # number of variables
s = 3 # nubmer of non-zero variables
X = matrix(rnorm(n*p), ncol=p)
beta = c(rep(3,s), rep(0,p-s))
y = 1 + X%*%beta + rnorm(n)
rm(n)
lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE)
lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE, penalty=list(lambda=0.1))
