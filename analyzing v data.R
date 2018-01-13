library(magrittr)
library(ggplot2)
library(progress)
library(quantreg)
library(viridis)

# Reading data ------
datasets <- vector(5, mode = "list")
for(i in 1:5) {
  dat <- read.csv(paste("data/data_1_", i, ".csv", sep = ""))
  dat$dataset <- paste("dataset", i, sep = "_")
  datasets[[i]] <- dat
}
dat <- do.call("rbind", datasets)

# Plotting Data -----------
ggplot(dat) +
  geom_point(aes(x = x, y = y)) +
  facet_wrap(~ dataset, scales = "free") +
  theme_bw()
# ggsave(filename = "tex/firstScatter.pdf", plot = last_plot(), units = "in",
#        height = 3, width = 6)

# Fitting Models -----
# mixture model
emfit <- hhlm(formula = y ~ dataset/x - 1,
            varFormula = ~ abs(x),
            data = dat,
            iterations = 100, delta = 10^-4)

# iteratively reweighted least squares
trimfit <- irwtls(y ~ x, data = dat,
                  varVariable = abs(dat$x),
                  dataset_identifier = dat$dataset,
                  iterations = 200, delta = 10^-4)

# quantile regression
quantfit <- rq(y ~ dataset/x - 1, data = dat)

# naive fit
naivefit <- lm(y ~ dataset/x - 1, data = dat)

# Coefficient estimates -----
emcoef <- coef(emfit$regModel) %>% matrix(ncol = 2) %>% t()
trimcoef <- sapply(trimfit$regModel, function(x) coef(x))
hybridcoef <- (trimcoef + emcoef) / 2
outtable <- t(hybridcoef)
colnames(outtable) <- c("b", "a")
outtable <- data.frame(outtable)
rownames(outtable) <- paste("data_1_", 1:5, ".csv", sep = "")
write.csv(outtable, "tex/vol_results.csv")

# Processing results ---------
trimFitted <- sapply(trimfit$regModel, function(x) x$fitted.values) %>% unlist()
estVarTrim <- trimfit$varEst
trimWeights <- sapply(trimfit$regModel, function(x) 1:length(x$residuals) %in% x$best) %>% unlist()
trimWeights <- trimWeights / (estVarTrim + 10^-4)

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
# Adding weights to graph
dataForPlot$weights <- with(dataForPlot, (1 - posterior) / (estVar + 10^-4))
dataForPlot$weights <- dataForPlot$weights / sum(dataForPlot$weights)
dataForPlot$trimWeights <- trimWeights
dataForPlot$trimWeights <- trimWeights / sum(trimWeights)
dataForPlot <- dataForPlot[order(abs(dataForPlot$x)), ]

# Plotting scatterplot
ggplot(dataForPlot) + geom_point(aes(x = x, y = y, col = trimWeights)) +
  geom_line(aes(x = x, y = yhat), col = "blue") +
  geom_line(aes(x = x, y = yhattrim), col = "red") +
  geom_line(aes(x = x, y = yhatQuant), col = "green") +
  geom_line(aes(x = x, y = yhatNaive), col = "orange") +
  facet_wrap(~ dataset, scales = "free") +
  scale_color_viridis() + theme_bw()
# ggsave(filename = "tex/emScatter.pdf", plot = last_plot(), units = "in",
#        height = 3, width = 6)

# Variance function plot ------
ggplot(dataForPlot) + geom_line(aes(x = (x), y = estVar), col = "blue") +
  geom_point(aes(x = (x), y = sqrt(estVarTrim)), col = "red") +
  geom_point(aes(x = x, y = sqrt(residuals^2))) + ylim(0, 40) + theme_bw() +
  ylab("squared-residuals") + xlab("x")
# ggsave(filename = "tex/varFunc.pdf", plot = last_plot(), units = "in",
#        height = 2.5, width = 5)

# Relative weights -----
plot(cumsum(sort(dataForPlot$weights, decreasing = TRUE)))

# Bootstrapping Estimate Variability ----------
# Some variables to help with saving results
ndatasets <- length(unique(dat$dataset))
datasetNum <- rep(1:ndatasets, 4)
method <- c("em", "trim", "quantile", "naive") %>%
  factor(., levels = .) %>%
  rep(5) %>% sort()

# Bootstrapping
set.seed(1)
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
  trimfit <- irwtls(y ~ x, data = bootDat,
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
bresults$method <- as.character(bresults$method)
bresults$method[bresults$method == "em"] <- "mixture"
ggplot(bresults) + geom_boxplot(aes(x = coef, y = estimate, col = method)) +
  facet_wrap(~ dataset, labeller = "label_both", scales = "free") +
  theme_bw()
# ggsave(filename = "tex/estBox.pdf", plot = last_plot(), units = "in",
#        height = 3, width = 6)

# Evaluating prediction error -----------
set.seed(1)
ndatasets <- length(unique(dat$dataset))
nCVs <- 20 # number of repeated CVs
nFolds <- 5 # number of folds in each repetition
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
    empred <- as.numeric((testX %*% coef(emfit$regModel)))
    emRMSE <- sqrt(mean((empred - test$y)^2))

    # trim
    trimfit <- irwtls(y ~ x, data = train,
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

    # mean prediction
    hybridPred <- (trimpred + empred) / 2
    hybridRMSE <- sqrt(mean((hybridPred - test$y)^2))

    # quantile
    quantfit <- rq(y ~ dataset/x - 1, data = train)
    quantRMSE <- sqrt(mean((predict(quantfit, newdata = test) - test$y)^2))

    # naive
    naivefit <- lm(y ~ dataset/x - 1, data = train)
    naiveRMSE <- sqrt(mean((predict(naivefit, newdata = test) - test$y)^2))

    # fold result
    foldresults[[f]] <- data.frame(iter = m,
                                   fold = f,
                                   method = c("em", "trim", "quantile", "naive", "hybrid"),
                                   rmse = c(emRMSE, trimRMSE, quantRMSE, naiveRMSE, hybridRMSE))
  }
  cvResults[[m]] <- do.call("rbind", foldresults)
}
close(pb)
cvResults <- do.call("rbind", cvResults)

# mean RMSE and SDs
library(dplyr)
library(xtable)
# Absolute RMSE
cvResults$method <- as.character(cvResults$method)
cvResults$method[cvResults$method == "em"] <- "mixture"
cvTable <- group_by(cvResults, method) %>%
  summarize(mRMSE = mean(rmse), rmseSD = sd(rmse)) %>%
  as.data.frame()

# Preparing latex table
rownames(cvTable) <- c("Hybrid", "Mixture", "Naive", "Quantile", "IRWTLS")
cvTable$method <- NULL
names(cvTable) <- c("RMSE", "SD")
xtable(t(cvTable), caption = "Mean RMSE in repeated cross-validation.")

# Errors Relative to Naive
relative <- subset(cvResults, method == "naive") %>%
  rename(naiveRMSE = rmse) %>% select(iter, fold, naiveRMSE) %>%
  merge(subset(cvResults, method != "naive"), all.x = TRUE, all.y = TRUE) %>%
  mutate(relError = log(rmse / naiveRMSE)) %>% group_by(method) %>%
  summarize(relRMSE = mean(relError), rmseSD = sd(relError)) %>%
  as.data.frame()
relative





