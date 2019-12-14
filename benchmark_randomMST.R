
#Install package
path <- getwd()
setwd(path)
#source("randomMstmv.R")
source("randomMstmv1.R")
#source("randomMstmv2.R")
install.packages("mstRforLR_1.2.zip", repos = NULL, type = "win.binary")
library(mstRforLR) 
install.packages("R.utils")
library(R.utils)
#################################################################################
#                                 Data call                                     #
#################################################################################
#Set pathway 
#path <- "C:/Users/monts/Box Sync/DIF in MST/Simulation MST/Data"
#setwd(path)


# Call true item parameters
itempar = read.csv("core.test.25.csv", header = T) #Test assemble/core.test.25.csv", header = T)
itembank = read.csv("item_bank_for_modules.csv", header = T)
Slope = itembank[, 4]
Location = itempar[, 3] # this is for selecting the DIF introduced to only 11%
it1 = data.frame(Slope, Location)
it1$c = 0
it1$d = 1
colnames(it1)[colnames(it1)=="Slope"] <- "a"
colnames(it1)[colnames(it1)=="Location"] <- "b"
it = it1

#################################################################################
#                              Assignment of modules                            #
#################################################################################
modules <- matrix(0, 72, 6)
modules[1:12, 1] <- modules[13:24, 2] <- modules[25:36, 3] <- modules[37:48, 4] <- modules[49:60, 5] <- modules[61:72, 6]<- 1

colSums(modules)

### Attach the module membership to each item
it_mod <- cbind(it, mod=0)
it_mod[1:12, "mod"]<- 1
it_mod[13:24, "mod"]<- 2
it_mod[25:36, "mod"]<- 3
it_mod[37:48, "mod"]<- 4
it_mod[49:60, "mod"]<- 5
it_mod[61:72, "mod"]<- 6

it_mod <- data.frame(it_mod)


# Cration of the transition matrix to define a 123 MST
trans <- matrix(0, 6, 6)
trans[1, 2:3] <- trans[2, 4:5] <- trans[3, 5:6]<-1
trans

plot.mst(trans)


# List for estimation of  theta

start<- list(fixModule = 1)
test1 <- list(moduleSelect="MFI", method = "EAP",  prob=c(1))
final<- list(method = "EAP")

#==================================================================
# First try to randomMST
#==================================================================


#Using microbenchmark in a smaller function
res1 <- randomMSTmv1(trueTheta = .99, itemBank = it, modules = modules, transMatrix = trans,
                    start = start, test = test1, final = final, allTheta = TRUE)
res1$selected.modules
res1$thetaProv
res1$seProv
res1$trueTheta
res1$thFinal

#==================================================================
# Second try to randomMST
#==================================================================


#################################################################################
#                                 MST SIMULATION                                #
#################################################################################

#===============================================================================
#                    ROUTING PROBABILITY = 1
#===============================================================================
## Time
path <- getwd() #"C:/Users/monts/Desktop/mstR/Results/"

Ntheta = 1000
NoReps <- c(10)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".25", 
                    DIFper = "11",
                    DIFloc = "Core", 
                    Prob = "1", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))




