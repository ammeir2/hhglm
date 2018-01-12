library(magrittr)

# Reading data ------
datasets <- vector(5, mode = "list")
for(i in 1:5) {
  dat <- read.csv(paste("data/data_1_", i, ".csv", sep = ""))
  dat$dataset <- paste("dataset", i, sep = "_")
  datasets[[i]] <- dat
}
dat <- do.call("rbind", datasets)

# Fitting Models -----
# mixture model
emfit <- hhlm(formula = y ~ dataset/x - 1,
            varFormula = ~ abs(x),
            data = dat,
            iterations = 100, delta = 10^-4)

# iteratively reweighted least squares
trimfit <- irtwls(y ~ x, data = dat,
                  varVariable = abs(dat$x),
                  dataset_identifier = dat$dataset,
                  iterations = 200, delta = 10^-4)

# quantile regression
quantfit <- rq(y ~ dataset/x - 1, data = dat)

# naive fit
naivefit <- lm(y ~ dataset/x - 1, data = dat)

# Processing results ---------
trimFitted <- sapply(trimfit$regModel, function(x) x$fitted.values) %>% unlist()
estVarTrim <- trimfit$varEst

# Plotting results ------------
dataForPlot <- data.frame(x = dat$x, y = dat$y,
                          yhat = emfit$predicted,
                          yhattrim = trimFitted,
                          yhatQuant = quantfit$fitted.values,
                          yhatNaive = predict(naivefit),
                          posterior = emfit$outlierPosteriors,
                          residuals = emfit$residuals,
                          dataset = dat$dataset,
                          estVar = sqrt(exp(predict(emfit$varModel))),
                          estVarTrim = estVarTrim)
dataForPlot$weights <- with(dataForPlot, (1 - posterior) / (estVar + 10^-3))
dataForPlot$weights <- dataForPlot$weights / sum(dataForPlot$weights)
dataForPlot$trimWeights <- trimWeights
dataForPlot$trimWeights <- trimWeights / sum(trimWeights)
dataForPlot <- dataForPlot[order(abs(dataForPlot$x)), ]

ggplot(dataForPlot) + geom_point(aes(x = x, y = y, col = weights)) +
  geom_line(aes(x = x, y = yhat), col = "blue") +
  geom_line(aes(x = x, y = yhattrim), col = "red") +
  geom_line(aes(x = x, y = yhatQuant), col = "green") +
  geom_line(aes(x = x, y = yhatNaive), col = "orange") +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_viridis() + theme_bw()
  # xlim(-10, 10) +
  # ylim(-15, 15)

ggplot(dataForPlot) + geom_line(aes(x = (x), y = estVar), col = "blue") +
  geom_line(aes(x = (x), y = sqrt(estVarTrim)), col = "red") +
  geom_point(aes(x = x, y = sqrt(residuals^2))) + ylim(0, 40)
plot(cumsum(sort(dataForPlot$weights, decreasing = TRUE)))

lm(y ~ dataset:x - 1, data = dat)

# Bootstrapping Estimate Variability ----------
# Some variables to help with saving results
ndatasets <- length(unique(dat$dataset))
datasetNum <- rep(1:ndatasets, 4)
method <- c("em", "trim", "quantile", "naive") %>%
  factor(., levels = .) %>%
  rep(5) %>% sort()

# doing bootstrap
library(progressBar)
bootReps <- 100
sepdat <- split(dat, dat$dataset)
pb <- txtProgressBar(min = 0, max = bootReps, style = 3)
bootResults <- vector(bootReps, mode = "list")
for(m in 1:bootReps) {
  setTxtProgressBar(pb, m)
  bootDat <- lapply(sepdat, function(x) x[sample.int(nrow(x), replace = TRUE), ]) %>%
    do.call("rbind", .)
  emfit <- hhlm(formula = y ~ dataset/x - 1,
                varFormula = ~ abs(x),
                data = bootDat,
                iterations = 100, delta = 10^-4)
  trimfit <- irtwls(y ~ x, data = bootDat,
                    varVariable = abs(bootDat$x),
                    dataset_identifier = bootDat$dataset,
                    iterations = 200, delta = 10^-4)
  trimcoef <- sapply(trimfit$regModel, coef)
  quantfit <- rq(y ~ dataset/x - 1, data = bootDat)
  naivefit <- lm(y ~ dataset/x - 1, data = bootDat)
  coefficients <- rbind(matrix(coef(emfit$regModel), ncol = 2),
                        t(trimcoef),
                        matrix(coef(quantfit), ncol = 2),
                        matrix(coef(naivefit), ncol = 2))
  iterResult <- data.frame(rep = m,
                           dataset = datasetNum,
                           method,
                           intercept = coefficients[, 1],
                           slope = coefficients[, 2])
  bootResults[[m]] <- iterResult
}
close(pb)

# Plotting bootstrap results
library(reshape2)
bresults <- do.call("rbind", bootResults)
bresults <- melt(bresults, id = c("rep", "dataset", "method"))
names(bresults)[4:5] <- c("coef", "estimate")
ggplot(bresults) + geom_boxplot(aes(x = coef, y = estimate, col = method)) +
  facet_wrap(~ dataset, labeller = "label_both", scales = "free")

# Evaluating prediction error -----------
ndatasets <- length(unique(dat$dataset))
nCVs <- 20
nFolds <- 5
cvResults <- vector(nCVs, mode = "list")
pb <- txtProgressBar(min = 0, max = nCVs * nFolds, style = 3)
pbInd <- 0
for(m in 1:nCVs) {
  folds <- split(dat, dat$dataset) %>%
    lapply(function(x) split(x, sample.int(nFolds, nrow(x), replace = TRUE)))
  foldresults <- vector(nFolds, mode = "list")
  for(f in 1:nFolds) {
    pbInd <- pbInd + 1
    setTxtProgressBar(pb, pbInd)

    # Getting test and train portions
    train <- folds %>% lapply(function(x) do.call("rbind", x[-f])) %>% do.call("rbind", .)
    test <- folds %>% lapply(function(x) do.call("rbind", x[f])) %>% do.call("rbind", .)

    # EM
    emfit <- hhlm(formula = y ~ dataset/x - 1,
                  varFormula = ~ abs(x),
                  data = train,
                  iterations = 100, delta = 10^-4)
    testX <- model.matrix(y ~ dataset/x - 1, data = test)
    emRMSE <- sqrt(mean((testX %*% coef(emfit$regModel) - test$y)^2))

    # trim
    trimfit <- irtwls(y ~ x, data = train,
                      varVariable = abs(train$x),
                      dataset_identifier = train$dataset,
                      iterations = 200, delta = 10^-4)
    trimpred <- vector(ndatasets, mode = "list")
    for(i in 1:ndatasets) {
      datasetInd <- test$dataset == paste("dataset", i, sep = "_")
      subX <- as.matrix(model.matrix(y ~ x, data = subset(test, datasetInd)))
      trimpred[[i]] <- subX %*% coef(trimfit$regModel[[i]])
    }
    trimpred <- unlist(trimpred)
    trimRMSE <- sqrt(mean((trimpred - test$y)^2))

    # quantile
    quantfit <- rq(y ~ dataset/x - 1, data = train)
    quantRMSE <- sqrt(mean((predict(quantfit, newdata = test) - test$y)^2))

    # naive
    naivefit <- lm(y ~ dataset/x - 1, data = train)
    naiveRMSE <- sqrt(mean((predict(naivefit, newdata = test) - test$y)^2))

    # fold result
    foldresults[[f]] <- data.frame(fold = f,
                                   method = c("em", "trim", "quantile", "naive"),
                                   rmse = c(emRMSE, trimRMSE, quantRMSE, naiveRMSE))
  }
  cvResults[[m]] <- do.call("rbind", foldresults)
}
close(pb)
cvResults <- do.call("rbind", cvResults)

# mean RMSE and SDs
library(dplyr)
library(xtable)
group_by(cvResults, method) %>%
  summarize(mRMSE = mean(rmse), rmseSD = sd(rmse)) %>%
  as.data.frame()




