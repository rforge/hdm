set.seed(1)
n = 100 #sample size
p = 100 # number of variables
s = 3 # nubmer of non-zero variables
X = matrix(rnorm(n*p), ncol=p)
beta = c(rep(3,s), rep(0,p-s))
y = 1 + X%*%beta + rnorm(n)

lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE)
lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE, penalty=list(lambda=0.1))

pred1 <- predict(lasso.reg)
pred2 <- predict(lasso.reg, newdata=X)
Xnew = matrix(rnorm(n/2*p), ncol=p)
pred1 <- predict(lasso.reg, newdata=Xnew)

data <- as.data.frame(cbind(y,X))
colnames(data) <- c("y", paste("V", 1:p, sep=""))
lasso.reg = rlasso(y~ V1+ V2 + V3 + V4 , post=TRUE, intercept=TRUE, data=data)
pred3 <- predict(lasso.reg)
pred4 <- predict(lasso.reg, newdata=data)

lasso.reg = rlasso(y~ V1+ V2 + V3 + V4 , post=TRUE, intercept=TRUE, data=data)
pred5 <- predict(lasso.reg)
pred6 <- predict(lasso.reg, newdata=data)

rhs <- paste("V", 1:100, sep="", collapse="+")
form <- as.formula(paste("y", "~", rhs))
lasso <- rlasso(form, data=data)
pred7 <- predict(lasso.reg)


#### Test significance test:

head(lasso.reg$dev)
lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE, model=FALSE)
summary(lasso.reg, all=FALSE)
lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE, model=TRUE)
summary(lasso.reg, all=FALSE)

lasso.reg = rlasso(y~X, post=TRUE, intercept=FALSE, model=TRUE)
summary(lasso.reg, all=FALSE)

set.seed(123456)
R <- 500
pvalues <- rep(NA, R)
n = 100 #sample size
p = 20 # number of variables
s = 0 # nubmer of non-zero variables
beta = c(rep(3,s), rep(0,p-s))
for (i in 1:R) {
X = matrix(rnorm(n*p), ncol=p)
y = 1 + X%*%beta + rnorm(n)
lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE, model=TRUE)
pvalues[i] <- summary(lasso.reg)$pval
}
