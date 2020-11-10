### runTMLE.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  6 2020 (11:32) 
## Version: 
## Last-Updated: Nov 10 2020 (08:46) 
##           By: Thomas Alexander Gerds
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
runTMLE <- function(no_cores,
                    K=10,#100              # end of follow-up (tau)
                    run.ltmle = FALSE,     # run ltmle, OR (!!)
                    run.ctmle = TRUE,      # run continuous tmle 
                    misspecify.init = FALSE,  # if TRUE, then misspecified model for outcome
                    M = 5, # 500           # number of simulation repetitions
                    n = 1000,              # sample size
                    max.iter=25,
                    seed,
                    progress.bar=3){

    #-------------------------------------------------------------------------------------------#
    ## repeat simulations (parallelize)
    #-------------------------------------------------------------------------------------------#
    message("\nEstimating psi with TMLE based on observed data:\n")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
        pb <- txtProgressBar(max = M,style = progress.bar,width=20)
    }
    out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
                   ) %dopar% {
                       if (progress.bar>0) { setTxtProgressBar(pb, m)}
                       repeat.fun(m,
                                  K=K,
                                  n=n,
                                  run.ltmle=run.ltmle,
                                  run.ctmle=run.ctmle,
                                  max.iter=max.iter,
                                  misspecify.init=misspecify.init,
                                  seed=seed+m)
                   }
    if (progress.bar>0) {    cat("\n")}
    if (no_cores>1) stopImplicitCluster()
    class(out) <- "watmle"
    out
}
######################################################################
### runTMLE.R ends here
