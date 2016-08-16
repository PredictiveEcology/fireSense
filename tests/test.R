library(SpaDES)

mySim <- simInit(
  times = list(start = 1, end = 2, timeunit = "year"),
  modules = list("fireSense"),
  paths = list(modulePath = " # replace with empty string instead"),
  inputs = data.frame(
    files = c("Z:/fireSenseProb.RData", "Z:/trajAgeMap.rds"),
    functions = c("load", "readRDS"),
    package = c("base", "base"),
    stringsAsFactors = FALSE)
)

spades(mySim)

hist(mySim$trajAgeMap$age)
