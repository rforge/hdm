HCHIV.lasso <- function(x,z,y,d,...) {
  n <- length(y)
  Z <- cbind(z,x)
  lasso.d.zx <- lasso(Z,d,...)
  lasso.y.x <- lasso(x,y,...)
  lasso.d.x <- lasso(x,d,...)
  if (sum(lasso.d.zx$index)==0) return(list(alpha=NA, se=NA))
  ind.dzx <- lasso.d.zx$index
  PZ <- Z[,ind.dzx]%*%solve(t(Z[,ind.dzx])%*%Z[,ind.dzx])%*%t(Z[,ind.dzx])%*%d
  lasso.PZ.x <- lasso(x,PZ,...)
  ind.PZx <- lasso.PZ.x$index
  Dr <- d- x[,ind.PZx]%*%solve(t(x[,ind.PZx])%*%x[,ind.PZx])%*%t(x[,ind.PZx])%*%PZ
  Yr <- lasso.y.x$residuals
  Zr <- lasso.PZ.x$residuals
  result <- tsls(Yr,Dr,x=NULL,Zr)
  return(list(coefficient=result$coefficient, se=sqrt(result$vcov)))
}

HCHIV.lasso.mult <- function(x,z,y,d,...) {
  #browser()
  d <- as.matrix(d)
  n <- dim(x)[1]
  d <- as.matrix(d)
  kd <- dim(d)[2]
  Z <- cbind(z,x)
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:kd, sep="")
  
  lasso.y.x <- lasso(x,y,...)
  Yr <- lasso.y.x$residuals
  Drhat <- NULL
  Zrhat <- NULL
  for (i in 1:kd) {
  lasso.d.x <- lasso(x,d[,i],...)
  lasso.d.zx <- lasso(Z,d[,i],...)
  if (sum(lasso.d.zx$index)==0) {
    Drhat <- cbind(Drhat, d[,i] - mean(d[,i]))
    Zrhat <- cbind(Zrhat, d[,i] - mean(d[,i]))
    next  
  }
  ind.dzx <- lasso.d.zx$index
  PZ <- Z[,ind.dzx,drop=FALSE]%*%solve(t(Z[,ind.dzx,drop=FALSE])%*%Z[,ind.dzx,drop=FALSE])%*%t(Z[,ind.dzx,drop=FALSE])%*%d[,i,drop=FALSE]
  lasso.PZ.x <- lasso(x,PZ,...)
  ind.PZx <- lasso.PZ.x$index
  Dr <- d[,i]- x[,ind.PZx,drop=FALSE]%*%solve(t(x[,ind.PZx,drop=FALSE])%*%x[,ind.PZx,drop=FALSE])%*%t(x[,ind.PZx,drop=FALSE])%*%PZ
  Zr <- lasso.PZ.x$residuals
  Drhat <- cbind(Drhat, Dr)
  Zrhat <- cbind(Zrhat, Zr)
  }
  result <- tsls(Yr,Drhat,x=NULL,Zrhat)
  coef <- as.vector(result$coefficient)
  se <- sqrt(diag(result$vcov))
  names(coef) <- names(se) <- colnames(d)
  return(list(coefficient= coef, se=se))
}
  
   