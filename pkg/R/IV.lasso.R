IV.lasso <- function(x,d,y,z, post=TRUE, ...) {
  
  d <- as.matrix(d)
  if(is.vector(x)) x <-  as.matrix(x)
  n <- length(y)
  kex <- dim(x)[2]
  ke <- dim(d)[2]
  kiv <- dim(z)[2]
  
  if (is.null(colnames(d))) colnames(d) <- paste("d", 1:ke, sep="")
  if (is.null(colnames(x)) & !is.null(x)) colnames(x) <- paste("x", 1:kex, sep="")
  # first stage regression
  Dhat <- NULL
  for (i in 1:ke) {
    di <- d[,i]
    lasso.fit <-lasso(z, di, post=post, ...)
    if (sum(lasso.fit$ind)==0) {
      dihat <- mean(di)
    } else {
      dihat <- z%*%lasso.fit$coefficients
    }
    Dhat <- cbind(Dhat, dihat)
  }
  Dhat <- cbind(Dhat, x)
  d <- cbind(d,x)
  # calculation coefficients
  alpha.hat <- solve(t(Dhat)%*%d)%*%(t(Dhat)%*%y)
  # calcualtion of the variance-covariance matrix
  residuals <- y - d%*%alpha.hat
  Omega.hat <- t(Dhat)%*%diag(as.vector(residuals^2))%*%Dhat #  Dhat.e <- Dhat*as.vector(residuals);  Omega.hat <- t(Dhat.e)%*%Dhat.e
  Q.hat.inv <- solve(t(d)%*%Dhat)
  vcov <- Q.hat.inv%*%Omega.hat%*%t(Q.hat.inv)
  rownames(alpha.hat) <- c(colnames(d))
  colnames(vcov) <- rownames(vcov) <- rownames(alpha.hat)
  return(list(coefficients=alpha.hat, vcov=vcov, residuals=residuals, samplesize=n))
}