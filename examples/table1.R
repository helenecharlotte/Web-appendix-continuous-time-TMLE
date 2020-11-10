### table1.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiesj Rijtgaard, Thomas Alexander Gerds
## Created: Nov  9 2020 (18:21) 
## Version: 
## Last-Updated: Nov 10 2020 (18:12) 
##           By: Thomas Alexander Gerds
##     Update #: 33
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

if (FALSE){ # be careful setting the number of cores and to watch the memory for large values of K   

    M <- 1000 # number of simulations
    nc <- 20 # number of cores
    n  <- 1000 # sample size

    # K=5
    table1.K5.true <- compute.true.psi0(n=100000,B=10,seed=17,K=5,verbose=2,no_cores=nc)
    saveRDS(table1.K5.true,file="./examples/table1-K5-true.rds")

    table1.K5.ltmle <- runLTMLE(no_cores=nc,K=5,misspecify.init = FALSE,M = M,n = n,seed=28,progress.bar=3)
    saveRDS(table1.K5.ltmle,file="./examples/table1-K5-ltmle.rds")
    
    table1.K5.ctmle <- runTMLE(no_cores=nc,K=5,misspecify.init = FALSE,M = M,n = n,max.iter=25,eps=0.001,seed=28,progress.bar=3)
    saveRDS(table1.K5.ctmle,file="./examples/table1-K5-ctmle.rds")

    # K=30
    table1.K30.true <- compute.true.psi0(n=100000,B=10,seed=17,K=30,verbose=2,no_cores=nc)
    saveRDS(table1.K30.true,file="./examples/table1-K30-true.rds")

    table1.K30.ltmle <- runLTMLE(no_cores=nc,K=30,misspecify.init = FALSE,M = M,n = n,seed=28,progress.bar=3)
    saveRDS(table1.K30.ltmle,file="./examples/table1-K30-ltmle.rds")
    
    table1.K30.ctmle <- runTMLE(no_cores=nc,K=30,misspecify.init = FALSE,M = M,n = n,max.iter=25,eps=0.001,seed=28,progress.bar=3)
    saveRDS(table1.K30.ctmle,file="./examples/table1-K30-ctmle.rds")
}

######################################################################
### table1.R ends here
