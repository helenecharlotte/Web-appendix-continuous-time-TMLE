### load.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  9 2020 (18:23) 
## Version: 
## Last-Updated: Nov 11 2020 (10:47) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
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
## source code
#-------------------------------------------------------------------------------------------#
source("./R/sim-data.R")
source("./R/compute-true.R")
source("./R/conTMLE.R")
source("./R/runTMLE.R")
source("./R/runLTMLE.R")
source("./R/helper-functions.R")
######################################################################
### load.R ends here
