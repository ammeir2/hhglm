# Reading data ------
datasets <- vector(5, mode = "list")
for(i in 1:5) {
  dat <- read.csv(paste("data/data_1_", i, ".csv", sep = ""))
  dat$dataset <- paste("dataset", i, sep = "_")
  datasets[[i]] <- dat
}
dat <- do.call("rbind", datasets)

fit <- hhlm(formula = y ~ dataset/x - 1,
            varFormula = ~ I(x^2),
            data = dat,
            iterations = 10, delta = 10^-2)

# plotting results -----
library(ggplot2)
library(viridis)

weights <- rep(1, nrow(dat))
for(m in 1:5) {
  trimfit <- ltsReg(y ~ dataset/x - 1, data = dat, weights = weights)
  residuals <- trimfit$residuals
  varmodel <- rq(log(residuals^2) ~ abs(x), data = dat)
  estSD <- sqrt(exp(varmodel$fitted.values))
  weights <- 1/(estSD^2 + 10^-4)
}

quantfit <- rq(y ~ dataset/x - 1, data = dat)
dataForPlot <- data.frame(x = dat$x, y = dat$y,
                          yhat = fit$predicted,
                          yhattrim = trimfit$fitted.values,
                          yhatQuant = quantfit$fitted.values,
                          posterior = fit$outlierPosteriors,
                          residuals = fit$residuals,
                          dataset = dat$dataset,
                          estVar = sqrt(exp(predict(fit$varModel))),
                          estVarTrim = sqrt(exp(varmodel$fitted.values)))
dataForPlot$weights <- with(dataForPlot, (1 - posterior) / (estVar + 10^-3))
dataForPlot$weights <- dataForPlot$weights / sum(dataForPlot$weights)
dataForPlot <- dataForPlot[order(abs(dataForPlot$x)), ]

ggplot(dataForPlot) + geom_point(aes(x = x, y = y, col = weights)) +
  geom_line(aes(x = x, y = yhat), col = "blue") +
  geom_line(aes(x = x, y = yhattrim), col = "red") +
  geom_line(aes(x = x, y = yhatQuant), col = "green") +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_viridis() + theme_bw()

ggplot(dataForPlot) + geom_line(aes(x = (x), y = estVar), col = "blue") +
  geom_line(aes(x = (x), y = sqrt(estVarTrim)), col = "red") +
  geom_point(aes(x = x, y = sqrt(residuals^2))) + ylim(0, 20)
plot(cumsum(sort(dataForPlot$weights, decreasing = TRUE)))

lm(y ~ dataset:x - 1, data = dat)

