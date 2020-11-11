### runTMLE.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  6 2020 (11:32) 
## Version: 
## Last-Updated: Nov 11 2020 (10:48) 
##           By: Thomas Alexander Gerds
##     Update #: 38
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
runTMLE <- function(no_cores,
                    K=10,                     # end of follow-up (tau)
                    misspecify.init = FALSE,  # if TRUE, then misspecified model for outcome
                    M = 5,                    # number of simulation repetitions
                    n = 1000,                 # sample size
                    max.iter=25,
                    eps = 0.0001,
                    seed,
                    progress.bar=3){

    message("\nEstimating psi with TMLE based on observed data:\n")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
        pb <- txtProgressBar(max = M,style = progress.bar,width=20)
    }
    out <- foreach(m=1:M, .errorhandling="pass") %dopar% {
        if (progress.bar>0) { setTxtProgressBar(pb, m)}
        # generate data
        dt <- sim.data(n,seed=seed+m,censoring=TRUE,K=K)
        # estimate 
        est.psi.A0 <- conTMLE(dt, 
                              targeting=2, 
                              smooth.initial=TRUE,
                              max.iter=max.iter,
                              eps=eps,
                              intervention.A0=function(L0, A0) logit(1*(A0==0)),
                              intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                              misspecify.init=misspecify.init)
        est.psi.A1 <- conTMLE(dt, 
                              targeting=2, 
                              smooth.initial=TRUE,
                              max.iter=max.iter,
                              eps=eps,
                              intervention.A0=function(L0, A0) logit(1*(A0==1)),
                              intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==1)),
                              misspecify.init=misspecify.init)
        names(est.psi.A0) <- paste0(names(est.psi.A1),".A0")
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
