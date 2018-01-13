#' Linear Regression for Data with Heteroskedastic Heavy-Tailed Data
#'
#' @importFrom quantreg rq
#' @export
hhlm <- function(formula, varFormula, data, iterations = 200, delta = 10^-4,
                 tol = 10^-3) {
  # Initializing -------
  initfit <- lm(formula, data = data)
  residuals <- initfit$residuals
  varFormula <- update(varFormula, log(residuals^2) ~ .) # Modifying variance equation
  data$residuals <- residuals
  lmvar <- lm(varFormula, data = data) # initial variance estimates
  estVar <- exp(predict(lmvar))
  estTdf <- findTdf(residuals / sqrt(estVar), rep(1, length(residuals))) # t-df estimation
  normalProb <- 0.5 # initializing probability of non-outlier
  tProb <- 1 - normalProb
  prevEst <- coef(initfit)

  # Generalized EM Algorithm
  for(m in 1:iterations) {
    # E Step ---------
    normDens <- dnorm(residuals, 0, sqrt(estVar), log = TRUE)
    tDens <- dt(residuals / sqrt(estVar), df = estTdf, log = TRUE)
    normPost <- 1 / (1 + exp(tDens + log(tProb) - normDens - log(normalProb)))
    tPost <- 1 - normPost # posterior probabilities of outlier

    # M-Step --------
    data$weights <- normPost / (estVar + delta) # Regression weights
    lmfit <- lm(formula, weights = weights, data = data) # Regression fit
    residuals <- lmfit$residuals
    data$residuals <- residuals
    data$normPost <- normPost
    lmvar <- rq(varFormula, data = data) # Estimating variance function
    estVar <- exp(predict(lmvar))
    estTdf <- findTdf(residuals / sqrt(estVar), tPost) # t-df estimation
    normalProb <- mean(normPost) # probability of non-outlier
    tProb <- 1 - normalProb

    # Stopping rule
    diff <- sum(abs(prevEst - coef(lmfit)) / abs(prevEst))
    prevEst <- coef(lmfit)
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

# A function for etimating t degrees of freedom via Maximum Likelihood -----
findTdf <- function(residuals, weights) {
  llim <- 2
  ulim <- 200
  max <- optimize(f = function(r) weighted.mean(dt(residuals, r, log = TRUE), weights),
                  lower = llim, upper = ulim, maximum = TRUE)$maximum
  return(max)
}
