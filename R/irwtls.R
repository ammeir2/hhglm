#' Fitting iteratively reweighted least trimmed squares models
#'
#' @importFrom robustbase ltsReg
#' @importFrom isotone gpava
#' @export
irwtls <- function(formula, data,
                   varVariable,
                   dataset_identifier,
                   iterations = 200,
                   delta = 10^-4,
                   tol = 10^-3) {
  weights <- rep(1, nrow(data))
  data <- data[order(dataset_identifier), ]
  for(m in 1:iterations) {
    data$weights <- weights
    # Fitting a diffrent regression model to each dataset
    trimFits <- by(data, dataset_identifier, function(d) ltsReg(formula, data = d, weights = weights))
    fitted <- sapply(trimFits, function(x) x$fitted.values) %>% unlist()
    residuals <- data$y - fitted
    isofit <- gpava(varVariable, log(residuals^2), p = 1) # isotonic regression
    estSD <- sqrt(exp(isofit$x)) # computing weights
    weights <- 1/(estSD^2 + 10^-4)

    # Stopping rule
    if(m > 1) {
      newcoef <- sapply(trimFits, coef)
      diff <- sum(abs(newcoef - prevcoef) / abs(prevcoef))
      if(diff < tol) break
    } else {
      prevcoef <- sapply(trimFits, coef)
    }
  }

  # Reporting
  result <- list(regModel = trimFits,
                 varEst = estSD^2,
                 residuals = residuals,
                 iterations = m)
  class(result) <- "irtwls"
  return(result)
}
