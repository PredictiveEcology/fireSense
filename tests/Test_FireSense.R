library(raster)
library(SpaDES)

modulePath <- "~/Documents/GitHub/McIntire-lab/modulesPrivate/"
paths <- list(modulePath = modulePath)

rasterOptions(maxmemory = 1e8)

# We start by fitting the spread model using the fireSense_SpreadFit module. In
# order to fit this model we need data about the shape of the fire size
# distribution. We can predict it using the fireSense_SizePredict module and fit
# the model to our study region using the fireSense_SizeFit module.

start <- 2000
end <- 2010

inputs <- rbind(
  data.frame(
    objectName = "dataFireSense_SizeFit",
    file = "C:/Z/Contrats/Eliot/Data/dataFireSense_SizeFit.rds",
    fun = "readRDS",
    package = "base",
    loadTime = start
  ),
  data.frame(
    objectName = "dataFireSense_SizePredict",
    file = "C:/Z/Contrats/Eliot/Data/dataFireSense_SizePredict_RASTER.rds",
    fun = "readRDS",
    package = "base",
    loadTime = start
  )
)

# Define module parameters
parameters <- list(
  fireSense_SizeFit = list(
    formula = list(beta = fireSize ~ hw + ot + dt + wt + MDC_07,
                   theta = fireSize ~ hw + ot + dt + wt + MDC_07),
    a = 1,
    itermax = 5000,
    trace = 100
  ),
  fireSense_SizePredict = list(
    intervalRunModule = 1
  )
)

sim <- simInit(
  times = list(start = start, end = end, timeunit = "year"),
  modules = list("fireSense_SizeFit",
                 "fireSense_SizePredict"),
  paths = paths,
  inputs = inputs,
  params = parameters
)

sim <- spades(sim)
datafireSense_SpreadFit <- sim$fireSense_SizePredicted
TP_Beta <- stack(lapply(datafireSense_SpreadFit, "[[", "beta"))
TP_Theta <- stack(lapply(datafireSense_SpreadFit, "[[", "theta"))

# Then we can fit frequency, escape and spread models using respectively the
# fireSense_FrequencyFit, fireSense_EscapeFit, fireSense_SpreadFit module and
# predict using their companion modules (fireSense_*Predict) and finally run
# fireSense.

# Define from where and how data will be loaded in the simList environment
inputs <- rbind(
  data.frame(
    objectName = "dataFireSense_FrequencyFit",
    file = "C:/Z/Contrats/Eliot/Data/dataFireSense_FrequencyFit.rds",
    fun = "readRDS",
    package = "base",
    loadTime = 1
  ),
  data.frame(
    objectName = "dataFireSense_EscapeFit",
    file = "C:/Z/Contrats/Eliot/Data/dataFireSense_EscapeFit.rds",
    fun = "readRDS",
    package = "base",
    loadTime = 1
  ),
  data.frame(
    objectName = "fireLoc_FireSense_SpreadFit",
    file = "C:/Z/Contrats/Eliot/Data/fireLoc_FireSense_SpreadFit.shp",
    fun = "shapefile",
    package = "raster",
    loadTime = 1
  ),
  data.frame(
    objectName = "dataFireSense_SpreadPredict",
    file = "C:/Z/Contrats/Eliot/Data/dataFireSense_SizePredict_RASTER.rds",
    fun = "readRDS",
    package = "base",
    loadTime = 1
  )
)

# Define module parameters
parameters <- list(
  fireSense_FrequencyFit = list(
    formula = nFires ~ hw:MDC_07 + cn:MDC_07 + ot:MDC_07 + dt:MDC_07 - 1,
    family = MASS::negative.binomial(theta = 1, link = "identity"),
    ub = list(coef = c(1, 1, 1, 1)),
    trace = 100,
    itermax = 100,
    nTrials = 100
  ),
  fireSense_FrequencyPredict = list(
    f = 100
  ),
  fireSense_EscapeFit = list(
    formula = cbind(escaped, nFires) ~ MDC_07 + hw + dt + ot + wt
  ),
  fireSense_SpreadFit = list(
    formula = ~ TP_Beta + TP_Theta,
    data = c("TP_Beta", "TP_Theta"),
    lower = c(.05, .1, .01, .3, 0.001, 0.001, 0.001),
    upper = c(.4, 10, .2, 4, .3, .3, .3),
    nCores = 1,
    itermax = 25
  ),
  fireSense = list(
    mapping = list(
      ignProb = "ignitionProb",
      escProb = "fireSense_EscapePredicted",
      sprProb = "fireSense_SpreadPredicted"
    )
  )
)

objects <- c("TP_Beta", "TP_Theta")

sim <- simInit(
  times = list(start = 1, end = 2, timeunit = "year"),
  modules = list("fireSense_FrequencyFit",
                 "fireSense_FrequencyPredict",
                 "fireSense_EscapeFit",
                 "fireSense_EscapePredict",
                 "fireSense_SpreadFit",
                 "fireSense_SpreadPredict",
                 "fireSense"),
  paths = paths,
  inputs = inputs,
  objects = objects,
  params = parameters
)

spades(sim)

hist(mySim$trajAgeMap$age)
