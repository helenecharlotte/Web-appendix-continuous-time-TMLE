### show.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Rytgaard
## Created: May, 2020 
#----------------------------------------------------------------------
## 
### Commentary:
## Look at output from TMLE estimation.
##  
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#-------------------------------------------------------------------------------------------#
## packages and generic functions
#-------------------------------------------------------------------------------------------#
library(data.table)
library(zoo)
library(stringr)
library(ltmle)
library(parallel)
library(foreach)
library(doParallel)

#-------------------------------------------------------------------------------------------#
## source relevant scripts
#-------------------------------------------------------------------------------------------#
setwd("~/research/SoftWare/Web-appendix-continuous-time-TMLE/")
source("./R/sim-data.R")
source("./R/est-fun.R")
source("./R/fit-density.R")
source("./R/runTMLE.R")
source("./R/helper-functions.R")
source("./R/repeat-fun.R")
source("./R/compute-true.R")

#-------------------------------------------------------------------------------------------#
## Example of the empirical studies
#-------------------------------------------------------------------------------------------#

# e0 <- runTMLE(no_cores=1,K=10,run.ltmle = FALSE,run.ctmle2 = TRUE,misspecify.Q = FALSE,only.A0 = FALSE,M = 1,n = 100,B = 1,N = 100,seed=4034249,progress.bar=4)

print(system.time(e1 <- runTMLE(no_cores=20,K=10,run.ltmle = FALSE,run.ctmle2 = TRUE,misspecify.Q = FALSE,only.A0 = FALSE,M = 200,n = 1000,B = 200,N = 1000,seed=4034249,progress.bar-1)))
saveRDS(e1,"~/research/SoftWare/Web-appendix-continuous-time-TMLE/examples/e1.rds")
e2 <- runTMLE(no_cores=20,K=10,run.ltmle = FALSE,run.ctmle2 = TRUE,misspecify.Q = FALSE,only.A0 = FALSE,M = 200,n = 1000,B = 200,N = 1000,seed=4034247,progress.bar=-1)
saveRDS(e2,"~/research/SoftWare/Web-appendix-continuous-time-TMLE/examples/e2.rds")

# misspecified Q
print(system.time(e1a <- runTMLE(no_cores=20,K=10,run.ltmle = FALSE,run.ctmle2 = TRUE,misspecify.Q = TRUE,only.A0 = FALSE,M = 200,n = 1000,B = 200,N = 1000,seed=4034249,progress.bar-1)))
saveRDS(e1a,"~/research/SoftWare/Web-appendix-continuous-time-TMLE/examples/e1a.rds")
e2a <- runTMLE(no_cores=20,K=10,run.ltmle = FALSE,run.ctmle2 = TRUE,misspecify.Q = TRUE,only.A0 = FALSE,M = 200,n = 1000,B = 200,N = 1000,seed=4034247,progress.bar=-1)
saveRDS(e2a,"~/research/SoftWare/Web-appendix-continuous-time-TMLE/examples/e2a.rds")
## larger K and requires a large memory
if (FALSE){
    e1 <- runTMLE(no_cores=1,
                  K=100,
                  run.ltmle = FALSE,
                  run.ctmle2 = TRUE,
                  misspecify.Q = FALSE,
                  only.A0 = FALSE,
                  M = 500,
                  n= 1000)
}

## e1
## plot(e1)




