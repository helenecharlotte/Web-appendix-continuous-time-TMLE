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
## set working directory
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    setwd("/home/ifsv/jhl781/research/phd/berkeley/continuousTMLE/")
} else {
    setwd("~/research/phd/berkeley/continuousTMLE/")
}

#-------------------------------------------------------------------------------------------#
## source relevant scripts
#-------------------------------------------------------------------------------------------#

source("./R/sim-data.R")
source("./R/est-fun.R")
source("./R/fit-density.R")
source("./R/helper-functions.R")
source("./R/repeat-fun.R")
source("./R/compute-true.R")

#-------------------------------------------------------------------------------------------#
## Example of the empirical studies
#-------------------------------------------------------------------------------------------#

e1 <- runTMLE(no_cores=1,
              K=10,
              run.ltmle = FALSE,
              run.ctmle2 = TRUE,
              misspecify.Q = FALSE,
              only.A0 = FALSE,
              M = 5,
              n= 1000)

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

e1
plot(e1)




