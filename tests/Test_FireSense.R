library(SpaDES)

modulePath <- "~/Documents/GitHub/McIntire-lab/modulesPrivate/"

# Define from where and how data will be loaded in the simList environment
inputs <- rbind(
  data.frame(
    objectName = "dataFireSense_FrequencyFit",
    file = "C:/Z/Contrats/Pessiere/DataInputs/dataFireSense_FrequencyFit_foudre.rds",
    fun = "readRDS",
    package = "base",
    loadTime = 1
  ),
 data.frame(
    objectName = "dataFireSense_EscapeFit",
    file = "C:/Z/Contrats/Pessiere/DataInputs/dataFireSense_EscapeFit_foudre.rds",
    fun = "readRDS",
    package = "base",
    loadTime = 1
  ),
 data.frame(
  objectName = "vegMap",
  file = "webDB",
  fun = "readRDS",
  package = "base",
  loadTime = 1
 ),
 data.frame(
   objectName = "ageMap",
   file = "C:/Z/trajAgeMap.rds",
   fun = "readRDS",
   package = "base",
   loadTime = 1
 )
)

# Define module parameters
parameters <- list(
  fireSense_SpreadFit = list(
    formula = ~ MDC_68 + cn + dt + ot + wt,
    lower = c(.2, .1, .01, .3, 0.001, 0.001),
    upper = c(.5, 10, .2, 4, .3, .3),
    trace = 5,
    nCores = 1,
    itermax = 5
  ),
  fireSense = list(
    mapping = list(
      ignProb = "ignitionProb",
      escProb = "fireSense_EscapePredicted",
      sprProb = "fireSense_SpreadPredicted"
    )
  )
)

paths <- list(modulePath = modulePath)

mySim <- simInit(
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
  params = parameters
)

spades(mySim)

hist(mySim$trajAgeMap$age)
