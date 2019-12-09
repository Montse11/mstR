# Testing if it produces the same output 

setwd("C:/Users/monts/Desktop/mstR/")
source("benchmark_randomMST.R")
source("randomMSTmv.R")

set.seed(1992)
res1 <- randomMSTmv(trueTheta = .99, itemBank = it, modules = modules, transMatrix = trans,
                  start = start, test = test1, final = final, allTheta = TRUE)
res1$selected.modules
res1$thetaProv
res1$seProv
res1$trueTheta
res1$thFinal
