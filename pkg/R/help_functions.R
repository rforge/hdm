format.perc <- function(probs, digits)
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
        "%")

hetero_se <- function(x,e,XpXinv) {
  k <- dim(x)[2]
  n <- dim(x)[1]
  
  V <- t((x*((e^2)%*%t(rep(1,k)))))%*%x
  
  vhetero <- n/(n-k)*XpXinv%*%V%*%t(XpXinv)
  se <- sqrt(diag(vhetero))
  
  return(list(se=se, vhetero=vhetero, V=V))
}

# function for calculation of the errors after choosing the five variables with the highest  correlation

init_values <- function(X,y, number=5, intercept=TRUE) {
  corr <- cor(y,X)
  kx <- dim(X)[2]
  index <- order(corr, decreasing = T)[1:min(number,kx)]
  coefficients <- rep(0,kx)
  if (intercept==TRUE) {
    reg <- lm(y~X[,index])
    coefficients[index] <- coef(reg)[-1]
  } else {
    reg <- lm(y~ -1 + X[,index])
    coefficients[index] <- coef(reg)
  }
  res <- list(residuals=reg$residuals, coefficients=coefficients)
  return(res)
}

