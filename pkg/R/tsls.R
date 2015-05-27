tsls <-
function(y,d,x,z) {
  n <- length(y)
  a1 <- dim(d)[2]
  a2 <- dim(x)[2]
  if (is.null(x)) {a2 <- 0}
  if (is.vector(x)) {a2 <- 1}
  if (is.vector(d)) {a1 <- 1}
  k <- a1 + a2
  X <- cbind(d,x)
  Z <- cbind(z,x)
  
  Mxz <- t(X)%*%Z
  Mzz <- solve(t(Z)%*%Z)
  M <- solve(Mxz%*%Mzz%*%t(Mxz))
  
  b <- M%*%Mxz%*%Mzz%*%(t(Z)%*%y)
  Dhat <- Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b2 <-  solve(t(Dhat)%*%X)%*%(t(Dhat)%*%y)
  e <- y - X%*%b
  VC1 <- as.numeric((t(e)%*%e/(n-k)))*M
 return(list(coefficients=b, vcov=VC1, residuals=e))
}
