### runTMLE.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  6 2020 (11:32) 
## Version: 
## Last-Updated: Nov 12 2020 (08:15) 
##           By: Thomas Alexander Gerds
##     Update #: 43
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
runTMLE <- function(K=10,                     
                    n = 1000,                 
                    misspecify.init = FALSE,  
                    seed,                     
                    M = 5,                    
                    no_cores=1,
                    max.iter=25,  # maximal number of iterations
                    eps = 0.0001, # convergence criterion
                    progress.bar=3, # show progress
                    verbose=FALSE, # get messages when running the code
                    ... # arguments passed to sim.data
                    ){ 
    message("\nEstimating psi with TMLE based on observed data:\n")
    if (M==1) progress.bar <-  -1
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
        pb <- txtProgressBar(max = M,style = progress.bar,width=20)
    }
    out <- foreach(m=1:M, .errorhandling="pass") %dopar% {
        if (progress.bar>0) { setTxtProgressBar(pb, m)}
        # generate data
        dt <- sim.data(n,seed=seed+m,censoring=TRUE,K=K,...)
        # estimate
        est.psi.A0 <- conTMLE(dt, 
                              targeting=2, 
                              smooth.initial=TRUE,
                              max.iter=max.iter,
                              eps=eps, dropin=FALSE,
                              verbose=verbose, 
                              intervention.A0=function(L0, A0) logit(1*(A0==0)),
                              intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                              misspecify.init=misspecify.init)
        est.psi.A1 <- conTMLE(dt, 
                              targeting=2, 
                              smooth.initial=TRUE,
                              max.iter=max.iter,
                              eps=eps, dropin=FALSE,
                              verbose=verbose, 
                              intervention.A0=function(L0, A0) logit(1*(A0==1)),
                              intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==1)),
                              misspecify.init=misspecify.init)
        names(est.psi.A0) <- paste0(names(est.psi.A0),".A0")
        names(est.psi.A1) <- paste0(names(est.psi.A1),".A1")
        return(c(est.psi.A0, est.psi.A1))
    }
    if (progress.bar>0) {    cat("\n")}
    if (no_cores>1) stopImplicitCluster()
    class(out) <- "watmle"
    out
}
######################################################################
### runTMLE.R ends here
