format.perc <- function(probs, digits) paste(format(100 * probs, trim = TRUE, 
                                                    scientific = FALSE, digits = digits), "%")

# function for calculation of the errors after choosing the five
# variables with the highest correlation

init_values <- function(X, y, number = 5, intercept = TRUE) {
  suppressWarnings(corr <- cor(y, X))
  kx <- dim(X)[2]
  index <- order(corr, decreasing = T)[1:min(number, kx)]
  coefficients <- rep(0, kx)
  if (intercept == TRUE) {
    reg <- lm(y ~ X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)[-1]
  } else {
    reg <- lm(y ~ -1 + X[, index, drop = FALSE])
    coefficients[index] <- coef(reg)
  }
  coefficients[is.na( coefficients)] <- 0
  res <- list(residuals = reg$residuals, coefficients = coefficients)
  return(res)
}

re.escape <- function(strings){
  vals <- c("\\\\", "\\[", "\\]", "\\(", "\\)", 
            "\\{", "\\}", "\\^", "\\$","\\*", 
            "\\+", "\\?", "\\.", "\\|")
  replace.vals <- paste0("\\\\", vals)
  for(i in seq_along(vals)){
    strings <- gsub(vals[i], replace.vals[i], strings)
  }
  strings
}
