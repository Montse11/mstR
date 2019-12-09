randomMSTmv = function (trueTheta = NULL, 
                        itemBank, 
                        modules, 
                        transMatrix, 
                        model = NULL, 
                        responses = NULL, 
                        genSeed = NULL, 
                        start = list(fixModule = NULL, seed = NULL, theta = 0, D = 1), 
                        test = list(method = "BM", priorDist = "norm", priorPar = c(0, 1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
                                    moduleSelect = "MFI", constantPatt = NULL, cutoff = NULL, prob = 1, seed.prob = NULL, score.range = "all"), 
                        final = list(method = "BM", priorDist = "norm", priorPar = c(0, 1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), alpha = 0.05), 
                        allTheta = FALSE, 
                        save.output = FALSE, 
                        output = c("path", "name", "csv")) {
  if (is.null(trueTheta) & is.null(responses)) 
    (break)("Either 'trueTheta' or 'responses' argument must be supplied", call. = FALSE)
  if (!testListMST(start, type = "start")$test) 
    stop(testListMST(start, type = "start")$message, call. = FALSE)
  if (!testListMST(test, type = "test")$test) 
    stop(testListMST(test, type = "test")$message, call. = FALSE)
  if (!testListMST(final, type = "final")$test) 
    stop(testListMST(final, type = "final")$message, call. = FALSE)
  if (!is.null(responses)) 
    assigned.responses <- TRUE
  else assigned.responses <- FALSE
  #InternalMST starts here
  internalMST <- function() {
    startList <- list(fixModule = start$fixModule, seed = start$seed, theta = 0, D = 1)
    startList$theta <- ifelse(is.null(start$theta), 0, start$theta)
    startList$D <- ifelse(is.null(start$D), 1, start$D)
    start <- startList
    testList <- list(method = NULL, priorDist = NULL, priorPar = c(0,1), 
                     range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
                     moduleSelect = "MFI", constantPatt = NULL, cutoff = NULL, 
                     prob = NULL, seed.prob = NULL, score.range = "all")
    testList$method <- ifelse(is.null(test$method), "BM", test$method)
    testList$priorDist <- ifelse(is.null(test$priorDist), "norm", test$priorDist)
    if (!is.null(test$priorPar)) {
      testList$priorPar[1] <- test$priorPar[1]
      testList$priorPar[2] <- test$priorPar[2]
    }
    if (!is.null(test$range)) {
      testList$range[1] <- test$range[1]
      testList$range[2] <- test$range[2]
    }
    testList$D <- ifelse(is.null(test$D), 1, test$D)
    if (!is.null(test$parInt)) {
      testList$parInt[1] <- test$parInt[1]
      testList$parInt[2] <- test$parInt[2]
      testList$parInt[3] <- test$parInt[3]
    }
    testList$moduleSelect <- ifelse(is.null(test$moduleSelect), 
                                    "MFI", test$moduleSelect)
    if (!is.null(test$constantPatt)) 
      testList$constantPatt <- test$constantPatt
    if (!is.null(test$cutoff)) 
      testList$cutoff <- test$cutoff
    if (is.null(test$prob)) 
      testList$prob <- 1
    else testList$prob <- test$prob
    if (!is.null(test$seed.prob)) 
      testList$seed.prob <- test$seed.prob
    if (!is.null(test$score.range)) 
      testList$score.range <- test$score.range
    test <- testList
    finalList <- list(method = NULL, priorDist = NULL, priorPar = c(0, 1), range = c(-4, 4), D = 1, parInt = c(-4, 4, 33), 
                      alpha = 0.05)
    finalList$method <- ifelse(is.null(final$method), "BM", 
                               final$method)
    finalList$priorDist <- ifelse(is.null(final$priorDist), 
                                  "norm", final$priorDist)
    if (is.null(final$priorPar) == FALSE) {
      finalList$priorPar[1] <- final$priorPar[1]
      finalList$priorPar[2] <- final$priorPar[2]
    }
    if (!is.null(final$range)) {
      finalList$range[1] <- final$range[1]
      finalList$range[2] <- final$range[2]
    }
    finalList$D <- ifelse(is.null(final$D), 1, final$D)
    if (!is.null(final$parInt)) {
      finalList$parInt[1] <- final$parInt[1]
      finalList$parInt[2] <- final$parInt[2]
      finalList$parInt[3] <- final$parInt[3]
    }
    finalList$alpha <- ifelse(is.null(final$alpha), 0.05, final$alpha)
    final <- finalList
    if (test$method == "score" & is.null(test$cutoff)) 
      (break)("'cutoff' argument of 'test' list must be supplied when 'method' is 'score'", 
              call. = FALSE)
    pr0 <- startModule(itemBank = itemBank, modules = modules, 
                       transMatrix = transMatrix, model = model, fixModule = start$fixModule, 
                       seed = start$seed, theta = start$theta, D = start$D)
    ITEMS <- pr0$items
    ITEMS.PER.MOD <- length(ITEMS)
    PAR <- rbind(pr0$par)
    MODULE <- pr0$module
    BEST.MOD <- NULL
    if (!is.null(responses)) 
      PATTERN <- responses[ITEMS]
    else PATTERN <- genPattern(trueTheta, PAR, model = model, 
                               D = test$D, seed = genSeed)
    if (test$method != "score") 
      TH <- thetaEst(PAR, PATTERN, model = model, D = test$D, 
                     method = test$method, priorDist = test$priorDist, 
                     priorPar = test$priorPar, range = test$range, 
                     parInt = test$parInt, current.th = start$theta, 
                     constantPatt = test$constantPatt, bRange = range(itemBank[, 
                                                                               2]))
    else TH <- sum(PATTERN)
    if (test$method != "score") 
      SETH <- semTheta(TH, PAR, x = PATTERN, model = model, 
                       D = test$D, method = test$method, priorDist = test$priorDist, 
                       priorPar = test$priorPar, parInt = test$parInt, 
                       constantPatt = test$constantPatt)
    else SETH <- NA
    thProv <- TH
    enter.mst <- TRUE
    if (sum(transMatrix[MODULE, ]) == 0) 
      enter.mst <- FALSE
    if (!enter.mst) {
      if (final$method != "score") {
        finalEst <- thetaEst(PAR, PATTERN, model = model, 
                             D = final$D, method = final$method, priorDist = final$priorDist, 
                             priorPar = final$priorPar, range = final$range, 
                             parInt = final$parInt)
        seFinal <- semTheta(finalEst, PAR, x = PATTERN, 
                            model = model, D = final$D, method = final$method, 
                            priorDist = final$priorDist, priorPar = final$priorPar, 
                            parInt = final$parInt)
        confIntFinal <- c(finalEst - qnorm(1 - final$alpha/2) * 
                            seFinal, finalEst + qnorm(1 - final$alpha/2) * 
                            seFinal)
      }
      else {
        finalEst <- sum(PATTERN)
        seFinal <- NA
        confIntFinal <- c(NA, NA)
      }
      RES <- list(trueTheta = trueTheta, selected.modules = MODULE, 
                  items.per.module = ITEMS.PER.MOD, transMatrix = transMatrix, 
                  model = model, testItems = ITEMS, itemPar = PAR, 
                  pattern = PATTERN, thetaProv = TH, seProv = SETH, 
                  thFinal = finalEst, seFinal = seFinal, ciFinal = confIntFinal, 
                  genSeed = genSeed, startFixModule = start$fixModule, 
                  startSeed = start$seed, startTheta = start$theta, 
                  startD = start$D, startThStart = pr0$thStart, 
                  startSelect = start$startSelect, provMethod = test$method, 
                  provDist = test$priorDist, provPar = test$priorPar, 
                  provRange = test$range, provD = test$D, moduleSelect = test$moduleSelect, 
                  constantPattern = test$constantPatt, cutoff = test$cutoff, 
                  prob = test$prob, seed.prob = test$seed.prob, 
                  score.range = test$score.range, best.module = NA, 
                  finalMethod = final$method, finalDist = final$priorDist, 
                  finalPar = final$priorPar, finalRange = final$range, 
                  finalD = final$D, finalAlpha = final$alpha, save.output = save.output, 
                  output = output, allTheta = NULL, assigned.responses = assigned.responses)
    }
    else {
      STAGE <- 0
      repeat {
        STAGE <- STAGE + 1
        if (length(test$prob) == 1) 
          PROB <- test$prob
        else PROB <- test$prob[STAGE]
        pr <- nextModule(itemBank, modules = modules, 
                         transMatrix = transMatrix, model = model, current.module = MODULE[length(MODULE)], 
                         theta = thProv, out = ITEMS, x = PATTERN, cutoff = test$cutoff, 
                         criterion = test$moduleSelect, parInt = test$parInt, 
                         priorDist = test$priorDist, priorPar = test$priorPar, 
                         D = test$D, range = test$range, prob = PROB, 
                         seed.prob = test$seed.prob)
        ITEMS <- c(ITEMS, pr$items)
        ITEMS.PER.MOD <- c(ITEMS.PER.MOD, length(pr$items))
        BEST.MOD <- c(BEST.MOD, pr$best.module)
        PAR <- rbind(PAR, pr$par)
        if (!is.null(responses)) 
          PATTERN <- c(PATTERN, responses[pr$items])
        else PATTERN <- c(PATTERN, genPattern(trueTheta, 
                                              pr$par, model = model, D = test$D, seed = genSeed))
        if (test$method != "score") 
          thProv <- thetaEst(PAR, PATTERN, model = model, 
                             D = test$D, method = test$method, priorDist = test$priorDist, 
                             priorPar = test$priorPar, range = test$range, 
                             parInt = test$parInt, current.th = TH[length(TH)], 
                             constantPatt = test$constantPatt, bRange = range(itemBank[, 
                                                                                       2]))
        else {
          if (test$score.range == "all") 
            thProv <- sum(PATTERN)
          else {
            nr <- length(pr$items)
            thProv <- sum(PATTERN[(length(PATTERN) - 
                                     nr + 1):(length(PATTERN))])
          }
        }
        TH <- c(TH, thProv)
        if (test$method != "score") 
          seProv <- semTheta(thProv, PAR, x = PATTERN, 
                             model = model, D = test$D, method = test$method, 
                             priorDist = test$priorDist, priorPar = test$priorPar, 
                             parInt = test$parInt, constantPatt = test$constantPatt)
        else seProv <- NA
        SETH <- c(SETH, seProv)
        MODULE <- c(MODULE, pr$module)
        if (sum(transMatrix[MODULE[length(MODULE)], ]) == 
            0) 
          break
      }
      if (final$method != "score") {
        finalEst <- thetaEst(PAR, PATTERN, model = model, 
                             D = final$D, method = final$method, priorDist = final$priorDist, 
                             priorPar = final$priorPar, range = final$range, 
                             parInt = final$parInt, current.th = TH[length(TH)], 
                             constantPatt = test$constantPatt, bRange = range(itemBank[, 
                                                                                       2]))
        seFinal <- semTheta(finalEst, PAR, x = PATTERN, 
                            model = model, D = final$D, method = final$method, 
                            priorDist = final$priorDist, priorPar = final$priorPar, 
                            parInt = final$parInt)
        confIntFinal <- c(finalEst - qnorm(1 - final$alpha/2) * 
                            seFinal, finalEst + qnorm(1 - final$alpha/2) * 
                            seFinal)
      }
      else {
        finalEst <- sum(PATTERN)
        seFinal <- NA
        confIntFinal <- c(NA, NA)
      }
      RES <- list(trueTheta = trueTheta, selected.modules = MODULE, 
                  items.per.module = ITEMS.PER.MOD, transMatrix = transMatrix, 
                  model = model, testItems = ITEMS, itemPar = PAR, 
                  pattern = PATTERN, thetaProv = TH, seProv = SETH, 
                  thFinal = finalEst, seFinal = seFinal, ciFinal = confIntFinal, 
                  genSeed = genSeed, startFixModule = start$fixModule, 
                  startSeed = start$seed, startTheta = start$theta, 
                  startD = start$D, provMethod = test$method, provDist = test$priorDist, 
                  provPar = test$priorPar, provRange = test$range, 
                  provD = test$D, moduleSelect = test$moduleSelect, 
                  constantPattern = test$constantPatt, cutoff = test$cutoff, 
                  prob = test$prob, seed.prob = test$seed.prob, 
                  score.range = test$score.range, best.module = BEST.MOD, 
                  finalMethod = final$method, finalDist = final$priorDist, 
                  finalPar = final$priorPar, finalRange = final$range, 
                  finalD = final$D, finalAlpha = final$alpha, save.output = save.output, 
                  output = output, allTheta = NULL, assigned.responses = assigned.responses)
    }
    if (allTheta) {
      prov.th <- prov.se <- NULL
      for (k in 1:nrow(RES$itemPar)) {
        if (test$method != "score") {
          prov.par <- rbind(RES$itemPar[1:k, ])
          prov.th[k] <- thetaEst(prov.par, RES$pattern[1:k], 
                                 model = model, D = test$D, method = test$method, 
                                 priorDist = test$priorDist, priorPar = test$priorPar, 
                                 range = test$range, parInt = test$parInt, 
                                 constantPatt = test$constantPatt, bRange = range(itemBank[, 
                                                                                           2]))
          prov.se[k] <- semTheta(prov.th[k], prov.par, 
                                 RES$pattern[1:k], model = model, D = test$D, 
                                 method = test$method, priorDist = test$priorDist, 
                                 priorPar = test$priorPar, parInt = test$parInt, 
                                 constantPatt = test$constantPatt)
        }
        else {
          prov.th[k] <- sum(RES$pattern[1:k])
          prov.se[k] <- NA
        }
      }
      RES$allTheta <- cbind(prov.th, prov.se)
      colnames(RES$allTheta) <- c("th", "se")
    }
    class(RES) <- "mst"
    return(RES)
  }
  resToReturn <- internalMST()
  if (save.output) {
    if (output[1] == "path") 
      wd <- paste(getwd(), "/", sep = "")
    else wd <- output[1]
    if (output[3] == "csv") 
      fileName <- paste(wd, output[2], ".csv", sep = "")
    else fileName <- paste(wd, output[2], ".txt", sep = "")
    capture.output(resToReturn, file = fileName)
  }
  return(resToReturn)
}
#<environment: namespace:mstRforLR>