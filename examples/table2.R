### table2.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiesj Rijtgaard, Thomas Alexander Gerds
## Created: Nov 10 2020 (17:58) 
## Version: 
## Last-Updated: Nov 11 2020 (17:21) 
##           By: Thomas Alexander Gerds
##     Update #: 11
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
source("./examples/load.R")

# be careful setting the number of cores and to watch the memory for large values of K
M <- 1000 # number of simulations
nc <- 40 # number of cores
n  <- 1000 # sample size


if (FALSE) { # K=5
  
    table2.K5.conTMLE <- runTMLE(no_cores=nc,K=5,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.0000001,seed=1920,progress.bar=-1)
    saveRDS(table2.K5.conTMLE,file="./examples/table2-K5-conTMLE.rds")

}

if (FALSE) { # K=30

    table2.K30.conTMLE <- runTMLE(no_cores=nc,K=30,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.0000001,seed=9170,progress.bar=-1)
    saveRDS(table2.K30.conTMLE,file="./examples/table2-K30-conTMLE.rds")

}


if (FALSE) { # K=50
    
    table2.K50.conTMLE <- runTMLE(no_cores=nc,K=50,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.0000001,seed=9170,progress.bar=-1)
    saveRDS(table2.K50.conTMLE,file="./examples/table2-K50-conTMLE.rds")

}

if (FALSE) { # K=100
    
    table2.K100.conTMLE <- runTMLE(no_cores=nc,K=100,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.0000001,seed=9170,progress.bar=-1)
    saveRDS(table2.K100.conTMLE,file="./examples/table2-K100-conTMLE.rds")

}



######################################################################
### table2.R ends here
