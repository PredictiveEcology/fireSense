library(raster)
library(SpaDES)

modulePath <- normalizePath("..")
paths <- list(modulePath = modulePath)

# rasterOptions(maxmemory = 1e8)

# We start by fitting the spread model using the fireSense_SpreadFit module. In
# order to fit this model we need data about the shape of the fire size
# distribution. We can predict it using the fireSense_SizePredict module and fit
# the model to our study region using the fireSense_SizeFit module.

start <- 2000
end <- 2010

inputs <- rbind(
  data.frame(
    objectName = "dataFireSense_SizeFit",
    file = normalizePath("../inputs/dataFireSense_SizeFit.rds"),
    fun = "readRDS",
    package = "base",
    loadTime = start
  ),
  data.frame(
    objectName = "dataFireSense_SizePredict",
    file = normalizePath("../inputs/dataFireSense_SizePredict_RASTER.rds"),
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
    initialRunTime =  2001,
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

TP_Beta <- stack(lapply(sim$fireSense_SizePredicted, "[[", "beta"))
TP_Theta <- stack(lapply(sim$fireSense_SizePredicted, "[[", "theta"))

# Then we can fit the spread model using the fireSense_SpreadFit module
inputs <- data.frame(
  objectName = "fireLoc_FireSense_SpreadFit",
  file = normalizePath("../inputs/fireLoc_FireSense_SpreadFit.shp"),
  fun = "shapefile",
  package = "raster",
  loadTime = 1
)

objects <- c("TP_Beta", "TP_Theta")

# Define module parameters
parameters <- list(
  fireSense_SpreadFit = list(
    formula = ~ TP_Beta + TP_Theta,
    data = c("TP_Beta", "TP_Theta"),
    lower = c(.05, .1, .01, .3, 0.001, 0.001, 0.001),
    upper = c(.1, .2, .2, 4, .3, .3, .3),
    nCores = 1,
    trace = 1,
    itermax = 5
  )
)

sim <- simInit(
  times = list(start = 1, end = 1, timeunit = "year"),
  modules = list("fireSense_SpreadFit"),
  paths = paths,
  inputs = inputs,
  objects = objects,
  params = parameters
)

sim <- spades(sim)
fireSense_SpreadFitted <- sim$fireSense_SpreadFitted


# Then we can fit frequency and escape models using respectively the 
# fireSense_FrequencyFit and fireSense_EscapeFit modules and predict using their
# companion modules (fireSense_*Predict). Then predict spread probabilities 
# using fireSense_SpreadPredict and finally run fireSense.


# Define from where and how data will be loaded in the simList environment
inputs <- rbind(
  data.frame(
    objectName = "dataFireSense_FrequencyFit",
    file = normalizePath("../inputs/dataFireSense_FrequencyFit.rds"),
    fun = "readRDS",
    package = "base",
    loadTime = 2000
  ),
  data.frame(
    objectName = "dataFireSense_FrequencyPredict",
    file = normalizePath("../inputs/dataFireSense_SizePredict_RASTER.rds"),
    fun = "readRDS",
    package = "base",
    loadTime = 2000
  ),
  data.frame(
    objectName = "dataFireSense_EscapeFit",
    file = normalizePath("../inputs/dataFireSense_EscapeFit.rds"),
    fun = "readRDS",
    package = "base",
    loadTime = 2000
  )
)

TP_Beta <- setNames(unstack(TP_Beta), nm = 2000:2010)
TP_Theta <- setNames(unstack(TP_Theta), nm = 2000:2010)
objects <- c("fireSense_SpreadFitted", "TP_Beta", "TP_Theta")

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
  fireSense_EscapePredict = list(
    data = "dataFireSense_FrequencyPredict"
  ),
  fireSense_SpreadPredict = list(
    data = c("dataFireSense_FrequencyPredict", "TP_Beta", "TP_Theta")
  ),
  fireSense = list(
    mapping = list(
      ignitionProb = "fireSense_FrequencyPredicted",
      escapeProb = "fireSense_EscapePredicted",
      spreadProb = "fireSense_SpreadPredicted"
    )
  )
)

sim <- simInit(
  times = list(start = start, end = end, timeunit = "year"),
  modules = list("fireSense_FrequencyFit",
                 "fireSense_FrequencyPredict",
                 "fireSense_EscapeFit",
                 "fireSense_EscapePredict",
                 "fireSense_SpreadPredict",
                 "fireSense"),
  paths = paths,
  inputs = inputs,
  objects = objects,
  params = parameters
)

sim <- spades(sim)


hist(mySim$trajAgeMap$age)
