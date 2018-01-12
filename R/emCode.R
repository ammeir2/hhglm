#' Linear Regression for Data with Heteroskedastic Heavy-Tailed Data
#'
#' @export
hhlm <- function(formula, varFormula, data, family, iterations = 20, delta = 10^-4) {
  # Initializing with robust regression -------
  robustFit <- glmrob(formula, data = data, family = "gaussian")
  residuals <- robustFit$residuals
  varFormula <- update(varFormula, log(residuals^2) ~ .)
  lmvar <- lm(varFormula, data = data)
  estVar <- exp(predict(lmvar))
  estTdf <- findTdf(residuals / sqrt(estVar), rep(1, length(residuals)))
  normalProb <- 0.5 # initializing probability of non-outlier
  tProb <- 1 - normalProb

  # Generalized EM Algorithm
  for(m in 1:iterations) {
    # E Step ---------
    normDens <- dnorm(residuals, 0, sqrt(estVar), log = TRUE)
    tDens <- dt(residuals / sqrt(estVar), df = estTdf, log = TRUE)
    normPost <- 1 / (1 + exp(tDens + log(tProb) - normDens - log(normalProb)))
    tPost <- 1 - normPost

    # M-Step --------
    lmfit <- lm(formula, weights = normPost / (estVar + delta))
    residuals <- lmfit$residuals
    lmvar <- lm(varFormula, weights = normPost / (estVar + delta))
    estVar <- exp(predict(lmvar))
    estTdf <- findTdf(residuals / sqrt(estVar), tPost)
    normalProb <- mean(normPost)
    tProb <- 1 - normalProb
  }

  # Reporting --------
  results <- list(regModel = lmfit,
                  varModel = estVar,
                  tdf = estTdf,
                  residuals = residuals,
                  outlierProbability = tProb)
  return(results)
}


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
  sds <- sapply(x, function(obs) sqrt(exp(sum(coef * c(obs, abs(obs), obs^2)))))
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


