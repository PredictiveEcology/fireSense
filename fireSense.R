# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense",
  description = "A fire module sensitive to weather and land cover changes",
  keywords = c("fire", "weather", "land cover"),
  authors = c(person("Jean", "Marchal", email="jean.d.marchal@gmail.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.0.0"),
  spatialExtent = raster::extent(c(-820340.251000005, -110340.251000005, 208055.595999999, 658055.595999999)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense.Rmd"),
  reqdPkgs = list("data.table", "dplyr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("startTime", "numeric", 1, NA, NA, desc = "Simulation time at which to initiate fire disturbance"),
    defineParameter("GCM", "character", "", NA, NA, desc = "Probability that a burning cell will continue burning for 1 iteration"),
    defineParameter("its", "numeric", Inf, NA, NA, desc = "Maximum number of iterations for the spread algorithm"),
    defineParameter("persistprob", "numeric", 0, 0, 1, desc = "Probability that a burning cell will continue burning for 1 iteration")
  ),
  inputObjects = data.frame(
    objectName = c("vegMap", "weatherVoxels", "ignitModel", "escapeModel", "sizeModel"),
    objectClass = c("data.frame", "data.frame", "list", "glm", "list"),
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = c("locis", "pc_escape", "fireSize"),
    objectClass = c("list", "numeric", "list"),
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

# load("C:/Users/jemar270/Desktop/Data/Climate/CMIP5/rcp85_mdc/CMIP5_data_bcc-bcc_csm1_1-r1i1p1.RData")
# load("Z:/PHD/Projets/Chap8/Data/escape_model.RData")
# load("Z:/PHD/Projets/Chap8/Data/Ignit_Model.RData")
# load("Z:/PHD/Projets/Chap8/Data/Size_Model.RData")


## Scale between 0 and 1 with possibly arbitrary min and max values
scale01 <- function(x, min = min(x), max = max(x)) (x - min) / (max - min)

doEvent.fireSense = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    wrapLoad <- function(x){
      x <- load(x, envir = environment())
      get(x)
    }
    
    envir(sim)$weatherVoxels <- as.data.table(wrapLoad(paste0("/home/jemar270/Chap8/Guillimin/Data/CMIP5_data_", params(sim)$fireSense$GCM,".RData")))
    envir(sim)$ignitModel <- wrapLoad("/home/jemar270/Chap8/Guillimin/Data/Ignit_Model.RData")
    envir(sim)$escapeModel <- wrapLoad("/home/jemar270/Chap8/Guillimin/Data/escape_model.RData")
    envir(sim)$sizeModel <- wrapLoad("/home/jemar270/Chap8/Guillimin/Data/Size_Model.RData")
    
    envir(sim)$ignitModel$mu <- gsub(envir(sim)$ignitModel$mu, pattern = "EE", replacement = "O") %>%
      gsub(pattern = "fw", replacement = "\"MDC_JUL\"") %>%
      gsub(pattern = "data\\[,", replacement = ".[[") %>%
      gsub(pattern = "\"]", replacement = "\"]]") %>%
      gsub(pattern = "param", replacement = "envir(sim)$ignitModel$par") %>%
      parse(text = .)
    
    envir(sim)$sizeModel <- envir(sim)$sizeModel[unlist(lapply(envir(sim)$sizeModel, function(x) any(x$LAMBDA == "MDC_JUN") & any(x$THETA == "MDC_MJ")))][[1]]
    envir(sim)$sizeModel$mu <-  parse(text = gsub(gsub(gsub(envir(sim)$sizeModel$mu, pattern = "o\\$p", replacement = "envir(sim)$sizeModel$par"), pattern = "L\\[,", replacement = ".[["), pattern = "\"]", replacement = "\"]]"))
    envir(sim)$sizeModel$theta <- parse(text = gsub(gsub(gsub(envir(sim)$sizeModel$theta, pattern = "o\\$p", replacement = "envir(sim)$sizeModel$par"), pattern = "L\\[,", replacement = ".[["), pattern = "\"]", replacement = "\"]]"))
    
# 
#    checkObject(sim, name = "landscapeTemplate")
#    checkObject(sim, name = "vegMap")
#    checkObject(sim, name = "weatherVoxels")
#     checkObject(sim, name = "ignitModel")
#     checkObject(sim, name = "escapeModel")
#     checkObject(sim, name = "sizeModel")

    sim <- sim$fireSenseInit(sim)
    sim <- scheduleEvent(sim, params(sim)$fireSense$startTime, "fireSense", "burn")
    
  } else if (eventType == "burn") {
    
    sim <- sim$fireSenseIgnition(sim)
    sim <- sim$fireSenseUpdateSpreadProb(sim)
    sim <- sim$fireSenseBurn(sim)
    
    sim <- scheduleEvent(sim, time(sim) + 1L, "fireSense", "burn")
    
  } else {
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}


fireSenseInit <- function(sim) {
  envir(sim)$locis <- vector("list", end(sim) - start(sim) + 1L)
  envir(sim)$pc_escape <- vector("numeric", end(sim) - start(sim) + 1L)
  envir(sim)$fireSize <- vector("list", end(sim) - start(sim) + 1L)
  
  invisible(sim)
}

fireSenseIgnition <- function(sim){
  ## Compute the mean of the negative binomial for the given set of weather and vegetation conditions
  envir(sim)$vegMap[, ':=' (HW = sum(HW), CN = sum(CN), DIST = sum(DIST), O = sum(O), WATER = sum(WATER)), by = px_id]
  envir(sim)$vegMap[, CELL_ID := (((ceiling(px_id / 710000)-1)*71) + ceiling((px_id - (ceiling(px_id / 7100)-1) * 7100)/100))]
  
  gc()
  
  ignit <- inner_join(envir(sim)$vegMap, dplyr::select(filter(envir(sim)$weatherVoxels, YEAR == time(sim)), CELL_ID, MDC_JUL), by = "CELL_ID") %>%
    with(., eval(envir(sim)$ignitModel$mu))
  
  ignit <- ignit / 10000 ## 100km2 to 1ha
  ignit[ignit < 0] <- 0
  
  ignited <- envir(sim)$vegMap[rbinom(n = length(ignit), size = 1, prob = ignit) == 1L, px_id]
  
  if(NROW(tmp <- envir(sim)$vegMap[px_id %in% ignited,])){

    escaped <- inner_join(tmp, dplyr::select(filter(envir(sim)$weatherVoxels, YEAR == time(sim)), CELL_ID, MDC_JUN), by = "CELL_ID") %>%
      predict(envir(sim)$escapeModel, newdata = ., type = "response") > runif(NROW(tmp))
    
    envir(sim)$locis[[time(sim) - start(sim) + 1L]] <- ignited[escaped]
    
    envir(sim)$pc_escape[time(sim) - start(sim) + 1L] <- sum(escaped) / length(ignited) * 100
  }
  
  invisible(sim)
}

fireSenseUpdateSpreadProb <- function(sim){
  if(!is.null(envir(sim)$locis[[time(sim) - start(sim) + 1L]])){
    load("/home/jemar270/Chap8/Guillimin/Data/DATA_agg_FS.RData", envir = environment())
    DATA <- DATA[DATA$CAUSE == "Lightning",]
    rge_MDC_JUN <- range(DATA$MDC_JUN)
    rge_MDC_MJ <- range(DATA$MDC_MJ)
    
    ## Compute LAMBDA
    LAMBDAS <- inner_join(envir(sim)$vegMap, dplyr::select(filter(sim$weatherVoxels, YEAR == time(sim)), CELL_ID, MDC_JUN), by = "CELL_ID") %>%
      mutate(MDC_JUN = pmin(1,pmax(0, sim$scale01(MDC_JUN, rge_MDC_JUN[1], rge_MDC_JUN[2])))) %>%
      with(.,eval(sim$sizeModel$mu)) %>%
      unlist(use.names = FALSE)
    LAMBDAS <- 1/LAMBDAS
    
    gc()
    
    ## Compute THETA
    THETAS <- inner_join(envir(sim)$vegMap, dplyr::select(filter(sim$weatherVoxels, YEAR == time(sim)), CELL_ID, MDC_MJ), by = "CELL_ID") %>%
      mutate(MDC_MJ = pmin(1,pmax(0,sim$scale01(MDC_MJ, rge_MDC_MJ[1], rge_MDC_MJ[2])))) %>%
      with(.,eval(sim$sizeModel$theta)) %>%
      unlist(use.names = FALSE)
    THETAS[THETAS < 1] <- 1
    THETAS <- log(THETAS)
    
    gc()
    
    ## From LAMBDA and THETA compute Pspread
    envir(sim)$spreadProb <- envir(sim)$landscapeTemplate
    v <- .116266 + .373580 / ((1 + (.016614 * LAMBDAS * THETAS)^(-5.993729))^3.295368)
    v[is.na(v)] <- 0
    envir(sim)$spreadProb[envir(sim)$vegMap[,px_id]] <- v
  }
  invisible(sim)
}

fireSenseBurn <- function(sim){
  if(!is.null(envir(sim)$locis[[time(sim) - start(sim) + 1L]])){
    fires <- SpaDES::spread(envir(sim)$spreadProb, loci=envir(sim)$locis[[time(sim) - start(sim) + 1L]],
                            spreadProb=envir(sim)$spreadProb, mask=NULL, persistence=params(sim)$fireSense$persistprob, directions=8L, 
                            iterations=params(sim)$fireSense$its, mapID=TRUE, returnIndices = TRUE) #, lowMemory = TRUE)
    envir(sim)$vegMap[px_id %in% fires[["indices"]], age := 0]
    envir(sim)$fireSize[[time(sim) - start(sim) + 1L]] <- tabulate(fires[["eventID"]])  
  }
  envir(sim)
  
  invisible(sim)
}


