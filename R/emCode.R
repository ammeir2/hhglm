#' Linear Regression for Data with Heteroskedastic Heavy-Tailed Data
#'
#' @importFrom quantreg rq
#' @export
hhlm <- function(formula, varFormula, data, iterations = 200, delta = 10^-4,
                 tol = 10^-3) {
  # Initializing with robust regression -------
  initfit <- lm(formula, data = data)
  residuals <- initfit$residuals
  varFormula <- update(varFormula, log(residuals^2) ~ .)
  data$residuals <- residuals
  lmvar <- lm(varFormula, data = data)
  estVar <- exp(predict(lmvar))
  estTdf <- findTdf(residuals / sqrt(estVar), rep(1, length(residuals)))
  normalProb <- 0.5 # initializing probability of non-outlier
  tProb <- 1 - normalProb
  prevEst <- coef(initfit)

  # Generalized EM Algorithm
  hackConst <- 1
  for(m in 1:iterations) {
    # print(m)
    # print(c(coef = coef(lmfit)))
    # print(c(varcoef = coef(lmvar)))
    # print(c(tdf = estTdf))
    # print(c(normProb = normalProb))
    # print(c(resids = quantile(residuals)))

    # E Step ---------
    normDens <- dnorm(residuals, 0, sqrt(estVar), log = TRUE)
    tDens <- dt(residuals / sqrt(hackConst * estVar), df = estTdf, log = TRUE)
    normPost <- 1 / (1 + exp(tDens + log(tProb) - normDens - log(normalProb)))
    tPost <- 1 - normPost

    # M-Step --------
    data$weights <- normPost / (estVar + delta)
    lmfit <- lm(formula, weights = weights, data = data)
    residuals <- lmfit$residuals
    data$residuals <- residuals
    data$normPost <- normPost
    #lmvar <- lm(varFormula, weights = normPost, data = data)
    lmvar <- rq(varFormula, data = data)
    estVar <- exp(predict(lmvar))
    estTdf <- findTdf(residuals / sqrt(hackConst * estVar), tPost)
    normalProb <- mean(normPost)
    tProb <- 1 - normalProb

    diff <- sum(abs(prevEst - coef(lmfit)) / abs(prevEst))
    prevEst <- coef(lmfit)
    # print(diff)
    if(diff < tol) break
  }

  # Reporting --------
  results <- list(regModel = lmfit,
                  varModel = lmvar,
                  tdf = estTdf,
                  predicted = predict(lmfit),
                  residuals = residuals,
                  outlierPosteriors = tPost,
                  outlierProbability = tProb,
                  iterations = m)
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
  llim <- 2
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


