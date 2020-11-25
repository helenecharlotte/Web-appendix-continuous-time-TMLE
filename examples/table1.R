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
setwd("/home/ifsv/jhl781/research/phd/berkeley/Web-appendix-continuous-time-TMLE/")
source("./examples/load.R")

if (TRUE){ # be careful setting the number of cores and to watch the memory for large values of K   

    M <- 1000 # number of simulations
    nc <- 40 # number of cores
    n  <- 1000 # sample size

    if (FALSE) {
        if (TRUE) {

            # K=5
            print(paste0("K=5, n=", n, ", M=", M, " nc=", nc))
            table1.K5.true <- compute.true(n=100000,B=50,seed=1113,K=5,no_cores=nc,progress.bar=-1)
            saveRDS(table1.K5.true,file="./examples/table1-K5-true.rds")
    
            table1.K5.conTMLE <- runTMLE(no_cores=nc,K=5,misspecify.init = FALSE,M = M,n = n,max.iter=25,eps=0.0000001,seed=20111,progress.bar=-1)
            saveRDS(table1.K5.conTMLE,file="./examples/table1-K5-conTMLE.rds")

            table2.K5.conTMLE <- runTMLE(no_cores=nc,K=5,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.0000001,seed=1920,progress.bar=-1)
            saveRDS(table2.K5.conTMLE,file="./examples/table2-K5-conTMLE.rds")
            
            table1.K5.ltmle <- runLTMLE(no_cores=nc,K=5,M = M,n = n,seed=2021,progress.bar=-1)
            saveRDS(table1.K5.ltmle,file="./examples/table1-K5-ltmle.rds")

        }
    
        if (FALSE) {
            # K=30
            nc <- 30
            print(paste0("K=30, n=", n, ", M=", M, " nc=", nc))
            table1.K30.true <- compute.true(n=100000,B=50,seed=17,K=30,no_cores=nc,progress.bar=-1)
            saveRDS(table1.K30.true,file="./examples/table1-K30-true.rds")

            table1.K30.conTMLE <- runTMLE(no_cores=nc,K=30,misspecify.init = FALSE,M = M,n = n,max.iter=25,eps=0.0000001,seed=1019,progress.bar=-1)
            saveRDS(table1.K30.conTMLE,file="./examples/table1-K30-conTMLE.rds")

            table2.K30.conTMLE <- runTMLE(no_cores=nc,K=30,misspecify.init = TRUE,M = M,n = n,max.iter=25,eps=0.0000001,seed=9170,progress.bar=-1)
            saveRDS(table2.K30.conTMLE,file="./examples/table2-K30-conTMLE.rds")

            if (TRUE) {
                table1.K30.ltmle <- runLTMLE(no_cores=nc,K=30,M = M,n = n,seed=5019,progress.bar=-1)
                saveRDS(table1.K30.ltmle,file="./examples/table1-K30-ltmle.rds")
            } 
        }

    }

    if (TRUE) {
        print(paste0("LTMLE: K=50, n=", n, ", M=", M, ", nc=", nc))
        table1.K50.ltmle <- runLTMLE(no_cores=nc,K=50,M = M,n = n,seed=5019,progress.bar=-1)
        saveRDS(table1.K50.ltmle,file="./examples/table1-K50-ltmle.rds")
    } 

    if (TRUE) {

        # K=100
        print(paste0("K=100, n=", n, ", M=", M, ", nc=", nc))
        nc <- 10
        table1.K100.true <- compute.true(n=100000,B=50,seed=1113,K=100,no_cores=nc,progress.bar=-1)
        saveRDS(table1.K100.true,file="./examples/table1-K100-true.rds")

        nc <- 40
        
        table1.K100.conTMLE <- runTMLE(no_cores=nc,K=100,misspecify.init=FALSE,
                                       M=M, n=n, max.iter=25,
                                       eps=0.0000001,seed=3019,
                                       progress.bar=-1)

        saveRDS(table1.K100.conTMLE,file="./examples/table1-K100-conTMLE.rds")

        table2.K100.conTMLE <- runTMLE(no_cores=nc,K=100,misspecify.init=TRUE,
                                      M=M, n=n, max.iter=25,
                                      eps=0.0000001,seed=5401,
                                      progress.bar=-1)

        saveRDS(table2.K100.conTMLE,file="./examples/table2-K100-conTMLE.rds")
    }

    if (FALSE) {

        # K=50
        print(paste0("K=50, n=", n, ", M=", M, " nc=", nc))
        table1.K50.true <- compute.true(n=100000,B=50,seed=1113,K=50,no_cores=nc,progress.bar=-1)
        saveRDS(table1.K50.true,file="./examples/table1-K50-true.rds")

        table1.K50.conTMLE <- runTMLE(no_cores=nc,K=50,misspecify.init=FALSE,
                                      M=M, n=n, max.iter=25,
                                      eps=0.0000001,seed=3019,
                                      progress.bar=-1)

        saveRDS(table1.K50.conTMLE,file="./examples/table1-K50-conTMLE.rds")

        table2.K50.conTMLE <- runTMLE(no_cores=nc,K=50,misspecify.init=TRUE,
                                      M=M, n=n, max.iter=25,
                                      eps=0.0000001,seed=5401,
                                      progress.bar=-1)

        saveRDS(table2.K50.conTMLE,file="./examples/table2-K50-conTMLE.rds")
    }
}

######################################################################
### table1.R ends here
