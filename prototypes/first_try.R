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
# loessVar <- loess(log(residuals^2) ~ x)
lmvar <- lm(log(residuals^2) ~ abs(x) + x^2)
plot(dat$x, sqrt(residuals^2)) # Checking that spline estimate is Reasonable
predvar <- exp(predict(lmvar)[order(x)])
lines(sort(x), sqrt(predvar)) # Clearly influenced by outliers
estTdf <- 3
normalProb <- 0.5 # initializing probability of non-outlier
tProb <- 1 - normalProb

# Starting EM -------
iterations <- 1
delta <- 10^-3
for(m in 1:iterations) {
  # E Step ---------
  estVar <- exp(predict(lmvar))
  normDens <- dnorm(residuals, 0, sqrt(estVar), log = TRUE)
  tDens <- dt(residuals / sqrt(estVar), df = estTdf, log = TRUE)
  normPost <- 1 / (1 + exp(tDens + log(tProb) - normDens - log(normalProb)))
  tPost <- 1 - normPost

  # M-Step --------
  lmfit <- lm(y ~ x, weights = normPost / (estVar + delta))
  residuals <- lmfit$residuals
  # loessVar <- loess(log(residuals^2) ~ abs(x), weights = normPost)
  # estVar <- exp(predict(loessVar))
  lmvar <- lm(log(residuals^2) ~ abs(x) + x^2, weights = normPost / (estVar + delta))
  estVar <- exp(predict(lmvar))
  estTdf <- findTdf(residuals / sqrt(estVar), tPost)
  normalProb <- mean(normPost)
  laplaceProb <- 1 - normalProb
}

library(ggplot2)
library(viridis)
weights <- normPost / (estVar + delta)
weights <- weights / sum(weights)
resultDat <- data.frame(x = x, y = y,
                        residsq = residuals^2, estVar = estVar,
                        normPost = normPost, weights = weights)
resultDat <- resultDat[order(resultDat$x), ]
ggplot(resultDat) + geom_point(aes(x = x, y = y, col = (normPost))) +
  geom_abline(intercept = lmfit$coefficients[1], slope = lmfit$coefficients[2]) +
  scale_color_viridis()


ggplot(resultDat) + geom_point(aes(x = x, y = sqrt(residsq), col = normPost)) +
  geom_line(aes(x = x, y = sqrt(estVar))) +
  scale_color_viridis()

plot(1:length(x), cumsum(sort(resultDat$weights, decreasing = TRUE)))

