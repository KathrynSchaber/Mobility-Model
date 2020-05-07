###
###
###
###
###


######################## ######################## ######################## ######################## ######################## 
### sensitivity runs ###
######################## ######################## ######################## ######################## ######################## 


##### move on days 1-3 ####
#########
rm(list = ls())

source('functions.R')
library(igraph, lib="myR")
library(snowfall, lib="myR")
library(Biodem, lib="myR")
library(popbio, lib="myR")
library(spatstat.data, lib="myR")
library(spatstat, lib="myR")
library(stats4)
#######
######################## parameters for all model runs ###########################
T = 200 # number of time steps to run
Nf = 600 # number of feeding stations
Nl = 600 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito extrinsic incubation period
tau = 1 # host intrinsic incubation period
rho_max = 5 # rho is host recovery 
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 8, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 8, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.75 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y

pL = points.poisson(Nl)
XLcoord = pL$x # larval habitat coordinate X
YLcoord = pL$y # larval habitat coordinate Y

betaParam1 = 21; betaParam2 = 11
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
if(length(which(alphaUneven > 0.85)) > 0){
  alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
}
alphaUnevenGradient = alphaUneven
alphaUnevenGradient[order(XLcoord ^ 2 + YLcoord ^ 2), ] = sort(alphaUneven, decreasing = TRUE)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))
hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts
HrossMacdonald = makeHrossMacdonald(Nf, Nh, Hh)
UrossMacdonaldRossMacdonald = makeUrossMacdonald(HrossMacdonald)
UrossMacdonaldUnevenBiting = makeUunevenBiting(HrossMacdonald)



############## H from SN (HM2), U Uneven, yes_changing move change SN/H NETOWRK METRICS ################
m <- Nf
House_counts <- as.integer(Hh) 

source('simulateEpidemic_sensitivity.R')
Nreps = 200

alpha = alphaUneven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

while(length(which(SmfNew > 350)) != 0){
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
  alpha = alphaUneven
  L = matrix(0, Nf, Nl)
  F = matrix(0, Nl, Nf)
  L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
  F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
  
  M_larvae = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = round(estimateLarvalEquilib())
  SmlNew = round(estimateMosquitoEquilibLarval())
  SmfNew = round(estimateMosquitoEquilibFeeding())
  print(SmfNew)
}

source('variables.R')
source('quick_H_approx.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

data_yes_changing_sensitivity = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 5)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
sfLibrary(igraph, lib="myR")
data_yes_changing_sensitivity = sfSapply(1 : Nreps, function(ii) simulateEpidemic_sensitivity(humanMeanR))
save(list = ls(), file = 'yes_changing_100_sensitivity.RData')
sfStop()


















##### yes move, 70% asymp ####
rm(list = ls())

source('functions.R')
library(stats4)
library(igraph)
library(snowfall)
######################## parameters for all model runs ###########################
T = 200 # number of time steps to run
Nf = 600 # number of feeding stations
Nl = 600 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito extrinsic incubation period
tau = 1 # host intrinsic incubation period
rho_max = 5 # rho is host recovery 
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 8, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 8, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.75 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y

pL = points.poisson(Nl)
XLcoord = pL$x # larval habitat coordinate X
YLcoord = pL$y # larval habitat coordinate Y

betaParam1 = 21; betaParam2 = 11
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
if(length(which(alphaUneven > 0.85)) > 0){
  alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
}
alphaUnevenGradient = alphaUneven
alphaUnevenGradient[order(XLcoord ^ 2 + YLcoord ^ 2), ] = sort(alphaUneven, decreasing = TRUE)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))
hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts
HrossMacdonald = makeHrossMacdonald(Nf, Nh, Hh)
UrossMacdonaldRossMacdonald = makeUrossMacdonald(HrossMacdonald)
UrossMacdonaldUnevenBiting = makeUunevenBiting(HrossMacdonald)



############## H from SN (HM2), U Uneven, yes_changing move change SN/H NETOWRK METRICS ################
m <- Nf
House_counts <- as.integer(Hh) 

source('simulateEpidemic_full_70_asymp.R')
Nreps = 200

alpha = alphaUneven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

while(length(which(SmfNew > 350)) != 0){
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
  alpha = alphaUneven
  L = matrix(0, Nf, Nl)
  F = matrix(0, Nl, Nf)
  L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
  F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
  
  M_larvae = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = round(estimateLarvalEquilib())
  SmlNew = round(estimateMosquitoEquilibLarval())
  SmfNew = round(estimateMosquitoEquilibFeeding())
  print(SmfNew)
}

source('variables.R')
source('quick_H_approx.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

data_yes_changing_70 = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 120)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
sfLibrary(igraph)
data_yes_changing_70 = sfSapply(1 : Nreps, function(ii) simulateEpidemic_full_70_asymp(humanMeanR))
save(list = ls(), file = 'yes_changing_70.RData')
sfStop()































######################## ######################## ######################## ######################## ######################## 
### yes move changes ###
######################## ######################## ######################## ######################## ######################## 


##### yes move, yes presymp ####
rm(list = ls())

source('functions.R')
library(stats4)
library(igraph)
library(snowfall)
######################## parameters for all model runs ###########################
T = 200 # number of time steps to run
Nf = 600 # number of feeding stations
Nl = 600 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito extrinsic incubation period
tau = 1 # host intrinsic incubation period
rho_max = 5 # rho is host recovery 
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 8, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 8, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.75 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y

pL = points.poisson(Nl)
XLcoord = pL$x # larval habitat coordinate X
YLcoord = pL$y # larval habitat coordinate Y

betaParam1 = 21; betaParam2 = 11
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
if(length(which(alphaUneven > 0.85)) > 0){
  alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
}
alphaUnevenGradient = alphaUneven
alphaUnevenGradient[order(XLcoord ^ 2 + YLcoord ^ 2), ] = sort(alphaUneven, decreasing = TRUE)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))
hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts
HrossMacdonald = makeHrossMacdonald(Nf, Nh, Hh)
UrossMacdonaldRossMacdonald = makeUrossMacdonald(HrossMacdonald)
UrossMacdonaldUnevenBiting = makeUunevenBiting(HrossMacdonald)


############## H from SN (HM2), U Uneven, yes_changing move change SN/H NETOWRK METRICS ################
m <- Nf
House_counts <- as.integer(Hh) 

source('simulateEpidemic_full.R')
Nreps = 200

alpha = alphaUneven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

while(length(which(SmfNew > 350)) != 0){
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
  alpha = alphaUneven
  L = matrix(0, Nf, Nl)
  F = matrix(0, Nl, Nf)
  L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
  F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
  
  M_larvae = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = round(estimateLarvalEquilib())
  SmlNew = round(estimateMosquitoEquilibLarval())
  SmfNew = round(estimateMosquitoEquilibFeeding())
  print(SmfNew)
}

source('variables.R')
source('quick_H_approx.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

data_yes_changing = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 120)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
sfLibrary(igraph)
data_yes_changing = sfSapply(1 : Nreps, function(ii) simulateEpidemic_full(humanMeanR))
save(list = ls(), file = 'yes_changing_100.RData')
sfStop()






























##### yes move, no pre ####
rm(list = ls())

source('functions.R')
library(stats4)
library(igraph)
library(snowfall)
######################## parameters for all model runs ###########################
T = 200 # number of time steps to run
Nf = 600 # number of feeding stations
Nl = 600 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito extrinsic incubation period
tau = 1 # host intrinsic incubation period
rho_max = 4 # rho is host recovery 
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 8, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 8, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.75 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y

pL = points.poisson(Nl)
XLcoord = pL$x # larval habitat coordinate X
YLcoord = pL$y # larval habitat coordinate Y

betaParam1 = 21; betaParam2 = 11
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
if(length(which(alphaUneven > 0.85)) > 0){
  alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
}
alphaUnevenGradient = alphaUneven
alphaUnevenGradient[order(XLcoord ^ 2 + YLcoord ^ 2), ] = sort(alphaUneven, decreasing = TRUE)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))
hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts
HrossMacdonald = makeHrossMacdonald(Nf, Nh, Hh)
UrossMacdonaldRossMacdonald = makeUrossMacdonald(HrossMacdonald)
UrossMacdonaldUnevenBiting = makeUunevenBiting(HrossMacdonald)



############## H from SN (HM2), U Uneven, yes_changing move change SN/H NETOWRK METRICS ################
m <- Nf
House_counts <- as.integer(Hh) 

source('simulateEpidemic_no_presymp.R')
Nreps = 200

alpha = alphaUneven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

while(length(which(SmfNew > 350)) != 0){
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
  alpha = alphaUneven
  L = matrix(0, Nf, Nl)
  F = matrix(0, Nl, Nf)
  L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
  F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
  
  M_larvae = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = round(estimateLarvalEquilib())
  SmlNew = round(estimateMosquitoEquilibLarval())
  SmfNew = round(estimateMosquitoEquilibFeeding())
  print(SmfNew)
}

source('variables.R')
source('quick_H_approx.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

data_yes_changing_no_presymp = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 120)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
sfLibrary(igraph)
data_yes_changing_no_presymp = sfSapply(1 : Nreps, function(ii) simulateEpidemic_no_presymp(humanMeanR))
save(list=ls(), file = 'yes_changing_no_presymp.RData')
sfStop()

















######################## ######################## ######################## ######################## ######################## 
### no move changes ###
######################## ######################## ######################## ######################## ######################## 


##### No move, presymp ####
rm(list = ls())

source('functions.R')
library(stats4)
library(igraph)
library(snowfall)
######################## parameters for all model runs ###########################
T = 200 # number of time steps to run
Nf = 600 # number of feeding stations
Nl = 600 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito extrinsic incubation period
tau = 1 # host intrinsic incubation period
rho_max = 5 # rho is host recovery 
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 8, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 8, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.75 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y

pL = points.poisson(Nl)
XLcoord = pL$x # larval habitat coordinate X
YLcoord = pL$y # larval habitat coordinate Y

betaParam1 = 21; betaParam2 = 11
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
if(length(which(alphaUneven > 0.85)) > 0){
  alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
}
alphaUnevenGradient = alphaUneven
alphaUnevenGradient[order(XLcoord ^ 2 + YLcoord ^ 2), ] = sort(alphaUneven, decreasing = TRUE)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))
hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts
HrossMacdonald = makeHrossMacdonald(Nf, Nh, Hh)
UrossMacdonaldRossMacdonald = makeUrossMacdonald(HrossMacdonald)
UrossMacdonaldUnevenBiting = makeUunevenBiting(HrossMacdonald)


############## H from SN (HM2), U Uneven, no move change SN/H NETOWRK METRICS ################
m <- Nf
House_counts <- as.integer(Hh) 

source('simulateEpidemic_no_move.R')
Nreps = 200

alpha = alphaUneven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

while(length(which(SmfNew > 350)) != 0){
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
  alpha = alphaUneven
  L = matrix(0, Nf, Nl)
  F = matrix(0, Nl, Nf)
  L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
  F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
  
  M_larvae = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = round(estimateLarvalEquilib())
  SmlNew = round(estimateMosquitoEquilibLarval())
  SmfNew = round(estimateMosquitoEquilibFeeding())
  print(SmfNew)
}

source('variables.R')
source('quick_H_approx.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

data_no_move = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 120)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
sfLibrary(igraph)
data_no_move = sfSapply(1 : Nreps, function(ii) simulateEpidemic_no_move(humanMeanR))
save(list=ls(), file = 'no_move.RData')
sfStop()
















##### No move, no presymp ####
rm(list = ls())

source('functions.R')
library(stats4)
library(igraph)
library(snowfall)
######################## parameters for all model runs ###########################
T = 200 # number of time steps to run
Nf = 600 # number of feeding stations
Nl = 600 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito extrinsic incubation period
tau = 1 # host intrinsic incubation period
rho_max = 5 # rho is host recovery 
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 8, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 8, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.75 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y

pL = points.poisson(Nl)
XLcoord = pL$x # larval habitat coordinate X
YLcoord = pL$y # larval habitat coordinate Y

betaParam1 = 21; betaParam2 = 11
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
if(length(which(alphaUneven > 0.85)) > 0){
  alphaEven = matrix(betaMean, nrow = Nl, ncol = 1)
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
}
alphaUnevenGradient = alphaUneven
alphaUnevenGradient[order(XLcoord ^ 2 + YLcoord ^ 2), ] = sort(alphaUneven, decreasing = TRUE)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))
hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts
HrossMacdonald = makeHrossMacdonald(Nf, Nh, Hh)
UrossMacdonaldRossMacdonald = makeUrossMacdonald(HrossMacdonald)
UrossMacdonaldUnevenBiting = makeUunevenBiting(HrossMacdonald)


############## H from SN (HM2), U Uneven, no move change SN/H NETOWRK METRICS ################
m <- Nf
House_counts <- as.integer(Hh) 

source('simulateEpidemic_no_move_no_presymp.R')
Nreps = 200

alpha = alphaUneven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())

while(length(which(SmfNew > 350)) != 0){
  alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = 1)
  alpha = alphaUneven
  L = matrix(0, Nf, Nl)
  F = matrix(0, Nl, Nf)
  L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
  F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
  
  M_larvae = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
  M_larvaeNew = round(estimateLarvalEquilib())
  SmlNew = round(estimateMosquitoEquilibLarval())
  SmfNew = round(estimateMosquitoEquilibFeeding())
  print(SmfNew)
}

source('variables.R')
source('quick_H_approx.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

data_no_move_no_presymp = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 120)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
sfLibrary(igraph)
data_no_move_no_presymp = sfSapply(1 : Nreps, function(ii) simulateEpidemic_no_move_no_presymp(humanMeanR))
save(list=ls(), file = 'no_move_no_presymp.RData')
sfStop()

######################## ######################## ######################## ######################## ########################
######################## ######################## ######################## ######################## ######################## 