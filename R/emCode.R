# Logit and expit -----
#' @export
logit <- function(x) log(x / (1 - x))
#' @export
expit <- function(x) 1 / (1 + exp(-x))

# Laplace Density Function ----
# zero mean is assumed
#' @importFrom robustbase glmrob
#' @export
dlaplace <- function(x, rate = 1, log = FALSE) {
  if(log) {
    return(dexp(abs(x), rate, log) + log(0.5))
  } else {
    return(0.5 * dexp(abs(x), rate, log))
  }
}

# A function for finding the MLE rate of a laplace random variable -----
#' @export
findLaplaceRateMLE <- function(residuals, weights) {
  ulim <- 1 / min(abs(residuals))
  llim <- 1 / max(abs(residuals))
  max <- optimize(f = function(r) weighted.mean(dlaplace(residuals, r, TRUE), weights),
                  lower = llim, upper = ulim, maximum = TRUE)$maximum
  return(max)
}

# A function for etimating t degrees of freedom -----
#' @export
findTdf <- function(residuals, weights) {
  llim <- 10^-3
  ulim <- 200
  max <- optimize(f = function(r) weighted.mean(dt(residuals, r, log = TRUE), weights),
                  lower = llim, upper = ulim, maximum = TRUE)$maximum
  return(max)
}


# A function for estimating a quadratic variance function -----
tAndNormLogLik <- function(coef, x, residuals, tdf) {
  sds <- sqrt(log(sum(coef * c(1, abs(x), x^2))))
  z <- residuals / sds
  dens <- -sum(weights * dnorm(z) + (1 - weights) * dt(z, df = tdf))
  print(c(dens, coef))
  return(dens)
}

#' @export
findVarFunction <- function(x, residuals, tdf, normalWeights) {
  minimum <- optim(par = c(1, 0, 0), fn = tAndNormLogLik,
                   x = x, residuals = residuals, tdf = tdf)
  return(minimum$par)
}


