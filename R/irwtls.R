#' @importFrom robustbase ltsReg
#' @importFrom isotone gpava
#' @export
irtwls <- function(formula, data,
                   varVariable,
                   dataset_identifier,
                   iterations = 200,
                   delta = 10^-4,
                   tol = 10^-3) {
  weights <- rep(1, nrow(data))
  data <- data[order(dataset_identifier), ]
  for(m in 1:iterations) {
    data$weights <- weights
    trimFits <- by(data, dataset_identifier, function(d) ltsReg(formula, data = d, weights = weights))
    fitted <- sapply(trimFits, function(x) x$fitted.values) %>% unlist()
    residuals <- data$y - fitted
    isofit <- gpava(varVariable, log(residuals^2), p = 0.5)
    estSD <- sqrt(exp(isofit$x))
    weights <- 1/(estSD^2 + 10^-4)

    if(m > 2) {
      newcoef <- sapply(trimFits, coef)
      diff <- sum(abs(newcoef - prevcoef) / abs(prevcoef))
      if(diff < tol) break
    } else {
      prevcoef <- sapply(trimFits, coef)
    }
  }

  result <- list(regModel = trimFits,
                 varEst = estSD^2,
                 residuals = residuals,
                 iterations = m)
  class(result) <- "irtwls"
  return(result)
}
