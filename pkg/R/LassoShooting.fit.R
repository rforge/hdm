LassoShooting.fit <-
function(x, y, lambda, control=list(maxIter=1000, optTol=10^(-5), zeroThreshold=10^(-6)), XX=NULL, Xy=NULL, beta.start=NULL) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(XX)) (XX <- crossprod(x))
  if(is.null(Xy)) (Xy <- crossprod(x,y))
  # Start from the LS solution for beta
  #beta <- solve(XX+diag(as.vector(lambda))%*%diag(1,p))%*%Xy
  if (is.null(beta.start)) {
  beta <- ginv(XX+diag(as.vector(lambda))%*%diag(1,p))%*%Xy
  beta[is.nan(beta)] <- 0
  } else {
  beta <- beta.start
  }
  # Start the log
  wp <- beta
  m <- 1
  XX2 <- XX*2
  Xy2 <- Xy*2
  
  while (m < control$maxIter) {
  beta_old <- beta    
  for (j in 1:p) {
    # Compute the Shoot and Update the variable
    S0 <- sum(XX2[j,]*beta) - XX2[j,j]*beta[j] - Xy2[j]
    if (sum(is.na(S0))>=1) {
        beta[j] <-0
        next
      }
    
    if (S0 > lambda[j]) beta[j] <- (lambda[j] - S0)/XX2[j,j]
    if (S0 < -1*lambda[j]) beta[j] <- (-1*lambda[j] - S0)/XX2[j,j]
    if (abs(S0) <= lambda[j]) beta[j] <- 0 
    }
  # Update the log
  wp <- cbind(wp,beta)
  # Check termination
  if (sum(abs(beta-beta_old)) < control$optTol) {break}
  m <- m+1
  }
  w <- beta
  w[abs(w)<control$zeroThreshold] <- 0
  return(list(coefficients=w, coef.list=wp, num.it=m))
}
