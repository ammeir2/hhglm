# loading data ------------
data_num <- 1
dat <- read.csv(paste("data/data_1_", data_num, ".csv", sep = ""))

# Plotting data ----------------
plot(dat$x, dat$y)

# fit robust regression ------
robustFit <- glmrob(y ~ x, data = dat, family = "gaussian")

# initial estimates -----------
y <- dat$y
x <- dat$x
residuals <- robustFit$residuals
loessVar <- loess(log(residuals^2) ~ x)
plot(dat$x, sqrt(residuals^2)) # Checking that spline estimate is Reasonable
predvar <- exp(predict(loessVar)[order(x)])
lines(sort(x), sqrt(predvar)) # Clearly influenced by outliers
laplaceRate <- 1 / max(abs(residuals)) # initializing from a high rate
normalProb <- 0.5 # initializing probability of non-outlier
laplaceProb <- 1 - normalProb

# Starting EM -------
iterations <- 10
delta <- 10^-3
for(m in 1:iterations) {
  # E Step ---------
  estVar <- exp(predict(loessVar))
  normDens <- dnorm(residuals, 0, sqrt(estVar), log = TRUE)
  laplaceDens <- dlaplace(residuals, rate = laplaceRate, log = TRUE)
  normPost <- 1 / (1 + exp(laplaceDens + log(laplaceProb) - normDens - log(normalProb)))
  lapPost <- 1 - normPost

  # M-Step --------
  lmfit <- lm(y ~ x, weights = normPost / (estVar + delta))
  residuals <- lmfit$residuals
  loessVar <- loess(log(residuals^2) ~ x, weights = normPost)
  laplaceRate <- findLaplaceRateMLE(residuals, lapPost)
  normalProb <- expit(predict(loess(logit(pmax(normPost, delta)) ~ x)))
  laplaceProb <- 1 - normalProb
}

library(ggplot2)
weights <- normPost / (estVar + delta)
weights <- weights / sum(weights)
resultDat <- data.frame(x = x, y = y,
                        residsq = residuals^2, estVar = estVar,
                        normPost = normPost, weights = weights)
resultDat <- resultDat[order(resultDat$x), ]
ggplot(resultDat) + geom_point(aes(x = x, y = y, col = weights)) +
  geom_abline(intercept = lmfit$coefficients[1], slope = lmfit$coefficients[2])

ggplot(resultDat) + geom_point(aes(x = x, y = sqrt(residsq), col = normPost)) +
  geom_line(aes(x = x, y = sqrt(estVar)))

plot(1:length(x), cumsum(sort(resultDat$weights, decreasing = TRUE)))

