### table1.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiesj Rijtgaard, Thomas Alexander Gerds
## Created: Nov  9 2020 (18:21) 
## Version: 
## Last-Updated: Nov 11 2020 (17:23) 
##           By: Thomas Alexander Gerds
##     Update #: 42
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
setwd("~/research/SoftWare/Web-appendix-continuous-time-TMLE/")
source("./examples/load.R")

if (TRUE){ # be careful setting the number of cores and to watch the memory for large values of K   

    M <- 1000 # number of simulations
    nc <- 40 # number of cores
    n  <- 1000 # sample size
    if (FALSE){
    # K=5
    table1.K5.true <- compute.true(n=100000,B=10,seed=1953,K=5,no_cores=nc,progress.bar=-1)
    saveRDS(table1.K5.true,file="./examples/table1-K5-true.rds")

    table1.K5.ltmle <- runLTMLE(no_cores=nc,K=5,misspecify.init = FALSE,M = M,n = n,seed=1986,progress.bar=-1)
    saveRDS(table1.K5.ltmle,file="./examples/table1-K5-ltmle.rds")
    
    table1.K5.conTMLE <- runTMLE(no_cores=nc,K=5,misspecify.init = FALSE,M = M,n = n,max.iter=25,eps=0.001,seed=2019,progress.bar=-1)
    saveRDS(table1.K5.conTMLE,file="./examples/table1-K5-conTMLE.rds")

    # K=30
    table1.K30.true <- compute.true(n=100000,B=10,seed=17,K=30,no_cores=nc,progress.bar=-1)
    saveRDS(table1.K30.true,file="./examples/table1-K30-true.rds")

    table1.K30.ltmle <- runLTMLE(no_cores=nc,K=30,misspecify.init = FALSE,M = M,n = n,seed=2052,progress.bar=-1)
    saveRDS(table1.K30.ltmle,file="./examples/table1-K30-ltmle.rds")
    }
    table1.K30.conTMLE <- runTMLE(no_cores=nc,K=30,misspecify.init = FALSE,M = M,n = n,max.iter=25,eps=0.001,seed=2085,progress.bar=-1)
    saveRDS(table1.K30.conTMLE,file="./examples/table1-K30-conTMLE.rds")
}

######################################################################
### table1.R ends here
