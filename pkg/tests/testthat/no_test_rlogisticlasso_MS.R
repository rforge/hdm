n <- 100
p <- 10
s <- 5
library(glmnet)
library(hdm2)
beta <- c(rep(2,s), rep(0,p-s))
X <- matrix(rnorm(n*p), ncol=p)
ystar <- X%*%beta + rnorm(n)
y <- as.integer(ystar>=0)

l <- 5
logistic.lassoA <- glmnet(X,y, family="binomial", lambda=l/(2*n))
logistic.lassoB <- rlogisticlasso(X,y, penalty=list(lambda.start= l), post=FALSE)
logistic.lassoC <- rlogisticlasso(X,y, penalty=list(lambda.start= l), post=TRUE)
logistic.lassoD <- rlogisticlasso(X,y, post=FALSE)
glmtest <- glm(y~X, family = binomial(link = "logit"))

logistic.lassoE <- rlogisticlasso(y~X, post=FALSE)
logistic.lassoF <- rlogisticlasso(y~X, post=TRUE)
logistic.lassoA <- glmnet(X,y, family="binomial")

summary(logistic.lassoC)
print(logistic.lassoB)
coef(logistic.lassoB)
yhat <- predict(logistic.lassoB)

d <- X[,1]
X <- X[,-1]

rlassologitEffect <- rlassologitEffectone(X,y,d)

rl <-  rlassologitEffect(X,y, index=c(1,2,7))
print(rl)
summary(rl)
confint(rl)


###########################
set.seed(2)
n <- 250
p <- 100
px <- 10

X <- matrix(rnorm(n*p), ncol=p)
colnames(X) <- paste("Var",1:100,sep="")
beta <- c(rep(1,px), rep(0,p-px))
P <- exp(X %*% beta)/(1+exp(X %*% beta))
y <- numeric(length=250)
for(i in 1:n){
  y[i] <- sample(x=c(1,0),size=1,prob=c(P[i],1-P[i]))
}
frame <- as.data.frame(cbind(y,X))
colnames(frame) <- c("Res",paste("V",1:100,sep=""))

d <- X[,1,drop=F]
X <- X[,-1,drop=F]

rlassologitEffect <- rlassologitEffectone(X,y,d)

rl <-  rlassologitEffect(X,y, index=c(1,2,7,20))
print(rl)
summary(rl)
confint(rl)
