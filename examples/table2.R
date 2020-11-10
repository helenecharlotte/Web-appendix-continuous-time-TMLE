### table2.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiesj Rijtgaard, Thomas Alexander Gerds
## Created: Nov 10 2020 (17:58) 
## Version: 
## Last-Updated: Nov 10 2020 (18:13) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
setwd("~/research/SoftWare/Web-appendix-continuous-time-TMLE/")
source("./examples/load.R")

if (FALSE){ # be careful setting the number of cores and to watch the memory for large values of K   

    M <- 1000  # number of simulations
    nc <- 20   # number of cores
    n  <- 1000 # sample size
    
    table2.K30.ctmle <- runTMLE(no_cores=nc,K=30,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.001,seed=28,progress.bar=3)
    saveRDS(table2.K30.ctmle,file="./examples/table2-K30-ctmle.rds")
}


######################################################################
### table2.R ends here