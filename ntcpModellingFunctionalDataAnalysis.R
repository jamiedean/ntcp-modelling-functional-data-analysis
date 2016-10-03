library(ggplot2)
library(fda.usc)
library(caret)
require(reshape2)
library(rms)
library(corrplot)
library(CalibrationCurves)
library(glmnet)

##################################################################################################################################
# User selections
##################################################################################################################################

# Select organ-at-risk to use: OM (oral mucosa), PM (pharyngeal mucosa)
oar <- 'OM'
# Select toxicity: dysphagia, mucositis
toxicityName <- 'mucositis'
# Select which model to use: pls, pca, mlr (multivariable logistic regression)
# Non-functional penalised logistic regression is included as an alternative dimensionality reduction approach to compare with the functional data analysis methods
dimRed <- 'pls'
# Select whether to use spatial dose metrics: TRUE (dose-length and dose-circumference histograms are used), FALSE (dose-volume histogram is used)
spatial <- FALSE
# Select number of bootstrap replicates
numberBootstraps = 2000

##################################################################################################################################
# Pre-processing
##################################################################################################################################

# Read in clinical (e.g. age, sex, chemotherapy) and toxicity data
dataAll <- read.csv(paste(toxicityName, 'FDA.csv', sep = ''), head = FALSE)

names(dataAll)[29] <- 'toxicity'
names(dataAll)[2] <- 'definitiveRT'
names(dataAll)[3] <- 'male'
names(dataAll)[4] <- 'age'
names(dataAll)[5] <- 'indChemo'
names(dataAll)[6] <- 'noConChemo'
names(dataAll)[7] <- 'cisplatin'
names(dataAll)[8] <- 'carboplatin'
names(dataAll)[9] <- 'cisCarbo'
names(dataAll)[11] <- 'hypopharynxLarynx'
names(dataAll)[12] <- 'oropharynxOralCavity'
names(dataAll)[13] <- 'nasopharynxNasalCavity'
names(dataAll)[14] <- 'unknownPrimary'
names(dataAll)[15] <- 'parotid'
names(dataAll)[16] <- 'V020'
names(dataAll)[17] <- 'V040'
names(dataAll)[18] <- 'V060'
names(dataAll)[19] <- 'V080'
names(dataAll)[20] <- 'V100'
names(dataAll)[21] <- 'V120'
names(dataAll)[22] <- 'V140'
names(dataAll)[23] <- 'V160'
names(dataAll)[24] <- 'V180'
names(dataAll)[25] <- 'V200'
names(dataAll)[26] <- 'V220'
names(dataAll)[27] <- 'V240'
names(dataAll)[28] <- 'V260'
names(dataAll)[10] <- 'independentValidation'

# Select only RMH and not WashU data
#data <- dataAll[dataAll[, 'independentValidation'] == 0,]
data <- dataAll

# Plot Spearman correlation matrix
colourPalette <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                                    "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))
corrplot(cor(data[,c(-1, -10)], method = 'spearman'), method = 'color', type = 'lower', tl.col = 'black', na.label = '-', col = colourPalette(200))

if (oar == 'PM' & spatial == TRUE) {
  dlhData <- read.csv(paste(oar, 'longExtFDA.csv', sep = ''), head = FALSE)
  dchData <- read.csv(paste(oar, 'circExtFDA.csv', sep = ''), head = FALSE)
  patientIDs.intersect = intersect(intersect(dlhData[, 1], dchData[, 1]), data[, 1])
  dlh <- fdata(subset(dlhData, V1 %in% patientIDs.intersect, select = -1))
  plot(dlh)
  dch <- fdata(subset(dchData, V1 %in% patientIDs.intersect, select = -1))
  plot(dch)
} else {
  dvhData <- read.csv(paste(oar, 'dvhForFDA.csv', sep = ''), head = FALSE)
  patientIDs.intersect = intersect(dvhData[, 1], data[, 1])
  dvh <- fdata(subset(dvhData, V1 %in% patientIDs.intersect, select = -1))
  plot(dvh)
}

data <- subset(data, V1 %in% patientIDs.intersect, select = -1)

# If there is an independent external validation dataset separate this out from the training data
xTrain = which(data$independentValidation == 0)
#print(xTrain)

if (oar == 'PM' & spatial == TRUE) {
  trainData = list('df' = data[xTrain,],
                   'dlh' = dlh[xTrain,], 'dlh.derivative' = dlh.derivative[xTrain,],
                   'dch' = dch[xTrain,], 'dch.derivative' = dch.derivative[xTrain,])
  testData = list('df' = data[-xTrain,],
                  'dlh' = dlh[-xTrain,], 'dlh.derivative' = dlh.derivative[-xTrain,],
                  'dch' = dch[-xTrain,], 'dch.derivative' = dch.derivative[-xTrain,])
} else {
  trainData = list('df' = data[xTrain,],
                   'dvh' = dvh[xTrain,])
  testData = list('df' = data[-xTrain,],
                  'dvh' = dvh[-xTrain,])
}

# Plot dose-volume histigrams for the training and independent external validation datasets
#plot(trainData$dvh)
#plot(testData$dvh)

##################################################################################################################################
# Model training
##################################################################################################################################

# Generate functional principal component or functional partial least squares basis functions
set.seed(0)
if (dimRed == 'pls' & spatial == FALSE) {
  # Find optimal partial least squares components and penalisation strength using Bayesian information criterion
  dvh.pls.cv = fregre.pls.cv(dvh[xTrain], data$toxicity[xTrain], kmax = 5, lambda = TRUE, P = c(0, 0, 1), criteria = 'SIC', norm = FALSE)
  print(dvh.pls.cv)
  print(dvh.pls.cv$pls.opt)
  print(dvh.pls.cv$lambda.opt)
  dvh.basis.pls = create.pls.basis(dvh[xTrain], data$toxicity[xTrain], l = dvh.pls.cv$pls.opt, norm = FALSE)
  basis.x = list('dvh' = dvh.basis.pls)
} else if (dimRed == 'pls' & spatial == TRUE) {
  # Find optimal partial least squares components and penalisation strength using Bayesian information criterion
  dlh.pls.cv = fregre.pls.cv(dlh[xTrain], data$toxicity[xTrain], kmax = 5, lambda = TRUE, P = c(0, 0, 1), criteria = 'SIC', norm = FALSE)
  dlh.basis.pls = create.pls.basis(dlh[xTrain], data$toxicity[xTrain], l = dlh.pls.cv$pls.opt, norm = FALSE)
  dch.pls.cv = fregre.pls.cv(dch[xTrain], data$toxicity[xTrain], kmax = 5, lambda = TRUE, P = c(0, 0, 1), criteria = 'SIC', norm = FALSE)
  dch.basis.pls = create.pls.basis(dch[xTrain], data$toxicity[xTrain], l = dch.pls.cv$pls.opt, norm = FALSE)
  basis.x = list('dlh' = dlh.basis.pls, 'dch' = dch.basis.pls)
} else if (dimRed == 'pca') {
  # Find optimal prinipal components and penalisation strength using Bayesian information criterion
  dvh.pca.cv = fregre.pc.cv(dvh[xTrain], data$toxicity[xTrain], kmax = 5, lambda = TRUE, P = c(0, 0, 1), criteria = 'SIC', norm = FALSE)
  print(dvh.pca.cv$pc.opt)
  print(dvh.pca.cv$lambda.opt)
  dvh.basis.pca = create.pc.basis(dvh[xTrain], l = dvh.pca.cv$pc.opt, norm = FALSE)
  basis.x = list('dvh' = dvh.basis.pca)
}

# Fit logistic regression model using functional partial least squares or functional principal components basis functions (for the dose-volume histogram data) and clinical covariates
if (dimRed == 'pls' | dimRed == 'pca') {
  model.final = fregre.glm(toxicity ~ dvh + cisplatin + carboplatin + cisCarbo +
                           indChemo + definitiveRT + male + age +
                           hypopharynxLarynx + nasopharynxNasalCavity + unknownPrimary + parotid,
                           data = trainData, family = binomial(), basis.x = basis.x)
  oddsRatios = exp(model.final$coefficients)
  summary.glm(model.final)
# Fit non-functional penalised (using LASSO) logistic regression model
} else if (dimRed == 'mlr') {
  covariates.mlr.train <- data.matrix(trainData$df[,!names(trainData$df) %in% c('toxicity', 'oropharynxOralCavity', 'noConChemo', 'independentValidation')])
  model.final = cv.glmnet(covariates.mlr.train, trainData$df$toxicity, family = 'binomial', alpha = 1, type.measure = 'auc')
  lambda.best = model.final$lambda.min
  oddsRatios = exp(coef(model.final))
}

print('Odds ratios:')
print(oddsRatios)

##################################################################################################################################
# Internal validation
##################################################################################################################################

# ncol has to be equal to the number of covariates plus the intercept
if (dimRed == 'pls' & spatial == FALSE) {
  num.covariates <- 12 + 5
} else if (dimRed == 'pls' & spatial == TRUE) {
  num.covariates <- 12 + 10
} else if (dimRed == 'pca') {
  num.covariates <- 12 + 5
} else if (dimRed == 'mlr') {
  num.covariates <- 25
}

bootmodel.beta = matrix(data = 0.0, nrow = numberBootstraps, ncol = num.covariates)
if (dimRed == 'pls' & spatial == FALSE) {
  boot.dvh.pls1 = matrix(data = NA, nrow = numberBootstraps, ncol = 270)
  boot.dvh.pls2 = matrix(data = NA, nrow = numberBootstraps, ncol = 270)
  boot.dvh.pls3 = matrix(data = NA, nrow = numberBootstraps, ncol = 270)
} else if (dimRed == 'pls' & spatial == TRUE) {
  boot.dlh.pls1 = matrix(data = NA, nrow = numberBootstraps, ncol = 14)
  boot.dlh.pls2 = matrix(data = NA, nrow = numberBootstraps, ncol = 14)
  boot.dlh.pls3 = matrix(data = NA, nrow = numberBootstraps, ncol = 14)
  boot.dch.pls1 = matrix(data = NA, nrow = numberBootstraps, ncol = 14)
  boot.dch.pls2 = matrix(data = NA, nrow = numberBootstraps, ncol = 14)
  boot.dch.pls3 = matrix(data = NA, nrow = numberBootstraps, ncol = 14)
} else if (dimRed == 'pca') {
  boot.dvh.pca1 = matrix(data = NA, nrow = numberBootstraps, ncol = 270)
  boot.dvh.pca2 = matrix(data = NA, nrow = numberBootstraps, ncol = 270)
  boot.dvh.pca3 = matrix(data = NA, nrow = numberBootstraps, ncol = 270)
}

optimism = 0
for (i in 1:numberBootstraps) {
  print('Bootstrap replicate:')
  print(i)
  set.seed(i)
  bootdata.df = trainData$df[sample(nrow(trainData$df), nrow(trainData$df), replace = TRUE), ]
  set.seed(i)
  bootdata.dvh = trainData$dvh[sample(nrow(trainData$dvh), nrow(trainData$dvh), replace = TRUE), ]
  # Bootstrap and plot FPLS or FPCA components
  if (dimRed == 'pls' & spatial == FALSE) {
    bootsample = list('df' = bootdata.df, 'dvh' = bootdata.dvh)
    boot.pls.cv = fregre.pls.cv(bootsample$dvh, bootsample$df$toxicity, kmax = 5, lambda = TRUE, P = c(0, 0, 1), criteria = 'SIC', norm = FALSE)
    print('Penalisation strength selected by Bayesian information criterion:')
    print(boot.pls.cv$lambda.opt)
    print('Partial least squares components selected by Bayesian information criterion:')
    print(boot.pls.cv$pls.opt)
    boot.dvh.pls <- fdata2pls(bootsample$dvh, bootsample$df$toxicity, ncomp = 5, lambda = boot.pls.cv$lambda.opt, P = c(0, 0, 1), norm = FALSE)
    boot.dvh.pls1[i,] = boot.dvh.pls$rotation[1]$data
    boot.dvh.basis.pls = create.pls.basis(bootsample$dvh, bootsample$df$toxicity, l = boot.pls.cv$pls.opt, lambda = boot.pls.cv$lambda.opt, P = c(0, 0, 1), norm = FALSE)
    boot.basis.x = list('dvh' = boot.dvh.basis.pls)
  } else if  (dimRed == 'pls' & spatial == TRUE) { 
    bootsample = list('df' = bootdata.df, 'dlh' = bootdata.dlh, 'dch' = bootdata.dch)
    boot.dlh.pls <- fdata2pls(bootsample$dlh, bootsample$df$toxicity, ncomp = length(dlh.pls.cv$pls.opt), norm = FALSE)
    boot.dlh.pls1[i,] = boot.dlh.pls$rotation[1]$data
    boot.dlh.basis.pls = create.pls.basis(bootsample$dlh, bootsample$df$toxicity, l = dlh.pls.cv$pls.opt, norm = FALSE)
    boot.dch.pls <- fdata2pls(bootsample$dch, bootsample$df$toxicity, ncomp = length(dch.pls.cv$pls.opt), norm = FALSE)
    boot.dch.pls1[i,] = boot.dch.pls$rotation[1]$data
    boot.dch.basis.pls = create.pls.basis(bootsample$dch, bootsample$df$toxicity, l = dch.pls.cv$pls.opt, norm = FALSE)
    boot.basis.x = list('dlh' = dlh.basis.pls, 'dch' = dch.basis.pls)
  } else if (dimRed == 'pca') {
    bootsample = list('df' = bootdata.df, 'dvh' = bootdata.dvh)
    boot.pca.cv = fregre.pc.cv(bootsample$dvh, bootsample$df$toxicity, kmax = 5, lambda = TRUE, P = c(0, 0, 1), criteria = 'SIC', norm = FALSE)
    print('Penalisation strength selected by Bayesian information criterion:')
    print(boot.pca.cv$lambda.opt)
    print('Principal components selected by Bayesian information criterion:')
    print(boot.pca.cv$pc.opt)
    boot.dvh.pca <- fdata2pc(bootsample$dvh, ncomp = 5, lambda = boot.pca.cv$lambda.opt, P = c(0, 0, 1), norm = FALSE)
    boot.dvh.pca1[i,] = boot.dvh.pca$rotation[1]$data
    boot.dvh.basis.pca = create.pc.basis(bootsample$dvh, l = boot.pca.cv$pc.opt, lambda = boot.pca.cv$lambda.opt, norm = FALSE)
    boot.basis.x = list('dvh' = boot.dvh.basis.pca)
  }
  # Bootstrap logistic regression models
  if ((dimRed == 'pls' | dimRed == 'pca') & spatial == FALSE) {
    bootmodel = fregre.glm(toxicity ~ dvh +
                           definitiveRT + male + age + indChemo + cisplatin + carboplatin + cisCarbo +
                           hypopharynxLarynx + nasopharynxNasalCavity + unknownPrimary + parotid,
                           data = bootsample, family = binomial(), basis.x = boot.basis.x)
  } else if ((dimRed == 'pls' | dimRed == 'pca') & spatial == TRUE) {
    bootmodel = fregre.glm(toxicity ~ dlh + dch +
                           definitiveRT + male + age + indChemo + cisplatin + carboplatin + cisCarbo +
                           hypopharynxLarynx + nasopharynxNasalCavity + unknownPrimary + parotid,
                           data = bootsample, family = binomial(), basis.x = boot.basis.x)
  } else if (dimRed == 'mlr') {
    covariates.mlr.train <- data.matrix(trainData$df[,!names(trainData$df) %in% c('toxicity', 'oropharynxOralCavity', 'noConChemo', 'independentValidation')])
    bootsample.covariates <- data.matrix(bootdata.df[,!names(bootdata.df) %in% c('toxicity', 'oropharynxOralCavity', 'noConChemo', 'independentValidation')])
    bootsample.toxicity = bootdata.df$toxicity
    bootmodel = cv.glmnet(bootsample.covariates, bootsample.toxicity, family = 'binomial', alpha = 1, type.measure = 'auc')
    boot.lambda.best = bootmodel$lambda.min
  }
  if (dimRed == 'pca') {
    #bootmodel.beta[i, 1:(12 + length(boot.pca.cv$pc.opt))] = bootmodel$coefficients
    bootmodel.beta[i, 1:12] = bootmodel$coefficients[1:12]
    bootmodel.beta[i, (12 + boot.pca.cv$pc.opt)] = bootmodel$coefficients[13:(12 + (length(boot.pca.cv$pc.opt)))]
  } else if (dimRed == 'pls') {
    bootmodel.beta[i, 1:(12 + length(boot.pls.cv$pls.opt))] = bootmodel$coefficients
  }
  if (dimRed != 'mlr') {
    boot.pred <- predict.fregre.glm(bootmodel, bootsample)
    boot.val = val.prob(boot.pred, bootsample$df$toxicity, pl = FALSE)
    boot.original.pred <- predict.fregre.glm(bootmodel, trainData)
    boot.original.val = val.prob(boot.original.pred, trainData$df$toxicity, pl = FALSE)
    optimism = optimism + (boot.val - boot.original.val)
  } else if (dimRed == 'mlr') {
    bootmodel.beta[i,] = data.matrix(coef(bootmodel))
    boot.pred <- predict(bootmodel, bootsample.covariates, s = boot.lambda.best, type = 'response')
    boot.val = val.prob(boot.pred, bootsample.toxicity, pl = FALSE)
    boot.original.pred <- predict(bootmodel, covariates.mlr.train, s = boot.lambda.best, type = 'response')
    boot.original.val = val.prob(boot.original.pred, trainData$df$toxicity, pl = FALSE)
    optimism = optimism + (boot.val - boot.original.val)
  }
}

optimism = optimism/numberBootstraps

# Plot first partial least squares component
if (dimRed == 'pls' & spatial == FALSE) {
  dvh.pls1.melt <- melt(boot.dvh.pls1)

  p = ggplot() + theme_bw() + geom_line(data = dvh.pls1.melt, aes(x = Var2, y = value, group = Var1), alpha = 0.025) +
    theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab('Fractional dose (cGy)') + ylab('Loadings') +
    scale_x_continuous(breaks = round(seq(0, 260, by = 20), 1))
  print(p)
# Plot first partial least squares components
} else if (dimRed == 'pls' & spatial == TRUE) {
  dlh.pls1.melt <- melt(boot.dlh.pls1)
  dch.pls1.melt <- melt(boot.dch.pls1)
  
  p = ggplot() + theme_bw() + geom_line(data = dlh.pls1.melt, aes(x = Var2, y = value, group = L1), alpha = 0.025) +
    theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab('Fractional dose (cGy)') + ylab('Loadings') +
    scale_x_discrete(labels = c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2',
                                '1.4', '1.6', '1.8', '2.0', '2.2', '2.4', '2.6'))
  print(p)
  p = ggplot() + theme_bw() + geom_line(data = dch.pls1.melt, aes(x = Var2, y = value, group = L1), alpha = 0.025) +
    theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab('Fractional dose (cGy)') + ylab('Loadings') +
    scale_x_discrete(labels = c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2',
                                '1.4', '1.6', '1.8', '2.0', '2.2', '2.4', '2.6'))
  print(p)
# Plot first principal component
} else if (dimRed == 'pca') {
  dvh.pca1.melt <- melt(boot.dvh.pca1)
    
  p = ggplot() + theme_bw() + geom_line(data = abs(dvh.pca1.melt), aes(x = Var2, y = value, group = Var1), alpha = 0.025) +
    theme(text = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab('Fractional dose (cGy)') + ylab('Loadings') +
    scale_x_continuous(breaks = round(seq(0, 260, by = 20), 1))
  print(p)
}

# Confidence intervals for odds ratios
x = matrix(data = 1.0, nrow = num.covariates, ncol = 1)
bootmodel.or.ci = matrix(data = 1.0, nrow = num.covariates, ncol = 2)
for (j in 1:num.covariates) {
  x[j] = j
  bootmodel.beta[, j]
  # 95 % confidence intervals
  or.ci = quantile(exp(bootmodel.beta[, j]), c(0.025, 0.975), na.rm = TRUE)
  bootmodel.or.ci[j,] = or.ci
}

oddsRatios.ci <- data.frame(bootmodel.or.ci)

if (dimRed == 'pls' & spatial == FALSE) {
  rownames(oddsRatios.ci) <- c('Intercept', 'definitiveRT', 'male', 'age', 'indChemo', 'cisplatin', 'carboplatin', 'cisCarbo',
                    'hypopharynxLarynx', 'nasopharynxNasalCavity', 'unknownPrimary', 'parotid',
                    'dvh.PLS1', 'dvh.PLS2', 'dvh.PLS3', 'dvh.PLS4', 'dvh.PLS5')
} else if (dimRed == 'pls' & spatial == TRUE) {
  rownames(oddsRatios.ci) <- c('Intercept', 'definitiveRT', 'male', 'age', 'indChemo', 'cisplatin', 'carboplatin', 'cisCarbo',
                    'hypopharynxLarynx', 'nasopharynxNasalCavity', 'unknownPrimary', 'parotid',
                    'dlh.PLS1', 'dch.PLS1')
} else if (dimRed == 'pca') {
  rownames(oddsRatios.ci) <- c('Intercept', 'definitiveRT', 'male', 'age', 'indChemo', 'cisplatin', 'carboplatin', 'cisCarbo',
                    'hypopharynxLarynx', 'nasopharynxNasalCavity', 'unknownPrimary', 'parotid',
                    'dvh.PC1', 'dvh.PC2', 'dvh.PC3', 'dvh.PC4', 'dvh.PC5')
} else if (dimRed == 'mlr') {
  rownames(oddsRatios.ci) <- c('Intercept', 'definitiveRT', 'male', 'age', 'indChemo', 'cisplatin', 'carboplatin', 'cisCarbo',
                    'hypopharynxLarynx', 'nasopharynxNasalCavity', 'unknownPrimary', 'parotid',
                    'V020', 'V040', 'V060', 'V080', 'V100', 'V120', 'V140', 'V160', 'V180', 'V200', 'V220', 'V240', 'V260')
}

# Plot odds ratios 95% confidence intervals
p = ggplot(oddsRatios.ci, aes(x = rownames(oddsRatios.ci), y = oddsRatios.ci$X1)) + xlab('Variable') + ylab('Odds ratio') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = oddsRatios.ci$X1, ymax = oddsRatios.ci$X2), width = 0.1) + 
  geom_hline(aes(yintercept = 1), colour = '#BB0000', linetype = 'dashed') +
  scale_y_log10()
print(p)


if (dimRed != 'mlr') {
  performance.apparent = val.prob(predict.fregre.glm(model.final, trainData), trainData$df$toxicity)
} else if (dimRed == 'mlr') {
  performance.apparent = val.prob(predict(model.final, covariates.mlr.train, s = lambda.best, type = 'response'), trainData$df$toxicity)
}

internalValidation = performance.apparent - optimism
print('Internal validation:')
print(internalValidation)

print('Odds ratios 95% confidence intervals:')
print(oddsRatios.ci)

##################################################################################################################################
# External validation
##################################################################################################################################

# External valiadtion dataset avaiable for dysphagia, but not mucositis
if (toxicityName == 'dysphagia') {
  numTestPatients <- length(data$toxicity[-xTrain])
  print(numTestPatients)
  if (dimRed != 'mlr') {
    pred <- predict.fregre.glm(model.final, testData)
  } else if (dimRed == 'mlr') {
    covariates.mlr.test <- data.matrix(testData$df[,!names(testData$df) %in% c('toxicity', 'oropharynxOralCavity', 'noConChemo', 'independentValidation')])
    pred <- predict(model.final, covariates.mlr.test, s = lambda.best, type = 'response')
  }
  y.pred <- ifelse(pred > 0.5, 1, 0)
  y.test <- data$toxicity[-xTrain]
  table(y.pred, y.test)

  # Model calibration curve and performance metrics
  externalValidation = val.prob(pred, y.test, smooth = FALSE, riskdist = FALSE, legendloc = FALSE, statloc = FALSE)
  val.prob.ci.2(pred, y.test, smooth = FALSE, logistic.cal = TRUE, riskdist = 'predicted', statloc = FALSE)

  # Bootstrap confidence intervals for external validation
  print(externalValidation)
  slope = externalValidation['Slope']
  intercept = externalValidation['Intercept']
  print(c(slope, intercept))
}

##################################################################################################################################
# Save model
##################################################################################################################################

# Delete data from model to be provided online
model.final$data = 0
save(model.final, file = paste('SavedModels/', toxicityName, oar, dimRed, '.rda', sep = ''))