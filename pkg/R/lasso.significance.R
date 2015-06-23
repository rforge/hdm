lasso.significance <- function(x, y, index = c(1:ncol(x)), alpha=0.05, method="bootstrap", B=500, ...) {
  k <- p1 <- length(index)
  n <- dim(x)[1]
  if (is.null(colnames(x))) 
    colnames(x) <- paste("V", 1:dim(x)[2], sep = "")
  table <- matrix(NA, ncol = 2, nrow = k)
  rownames(table) <- colnames(x)[index]
  colnames(table) <- c("coeff.", "se.")
  reside <- matrix(NA, nrow=n, ncol=p1)
  residv <- matrix(NA, nrow=n, ncol=p1)
  for (i in 1:k) {
    d <- x[, index[i], drop = FALSE]
    Xt <- x[, -index[i], drop = FALSE]
    result <- try(col <- HC.lasso(Xt, y, d, I3 = NULL, ...))
    if (class(result) == "try-error") {
      next
    }
    else {
      table[i, 1] <- col$alpha
      table[i, 2] <- col$se
      reside[,i] <- col$residuals$epsilon
      residv[,i] <- col$residuals$v
    }
  }
  phi <- reside*residv
  m <- 1/sqrt(colMeans(phi^2))
  phi <- t(t(phi)/m)
  sigma <- sqrt(colMeans(phi^2))
  sim <- vector("numeric", length=B)
  for (i in 1:B) {
  xi <- rnorm(n)
  phi_temp <- phi*xi
  Nstar <- 1/sqrt(n)*colSums(phi_temp)
  sim[i] <- max(abs(Nstar))
  }
  hatc <- quantile(sim, probs=1-alpha)
  res <- matrix(NA, ncol = 4, nrow = k)
  rownames(res) <- colnames(x)[index]
  colnames(res) <- c("coeff.", "se", "lower CI", "upper CI")
  res[,1] <- table[,1]
  res[,2] <- table[,2]
  res[,3] <- res[,1] - hatc*1/sqrt(n)*sigma
  res[,4] <- res[,1] + hatc*1/sqrt(n)*sigma
  return(res)
}