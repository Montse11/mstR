# Measuring performance
path <- getwd()
setwd(paste0(path, "/mstR"))
#Call packages
#install.packages("profvis")
library(profvis)
source("benchmark_randomMST.R")

profvis({
  for (i in 1:1){ 
    set.seed(i*1992)
    th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                          rnorm(Ntheta, mean= -1.28, sd = .87),
                          rnorm(Ntheta, mean= -1.32, sd = .86)) 
    ressim <- NULL
    for (k in 1:NoReps[1]){
      test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(1), seed.prob = (i*k+47405))
      
      res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                        start = start, 
                        test = test2, 
                        final = final)
      
      data[k, c("rep")] <- i
      data[k, c("TH0")] <- res2$trueTheta
      data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
      data[k, c("OP1", "OP2")] <- res2$best.module
      data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
      data[k, c("th1", "th2", "th3")] <- res2$thetaProv
      data[k, c("seT")] <- res2$seFinal
      data[k, c("set1", "set2", "set3")] <- res2$seProv
      #ressim = cbind(res, data)
    }
    saveObject(data, file = paste0(path, "Core.25.11.p1", "_i", i, ".Rbin"))
    
  }
}
)

