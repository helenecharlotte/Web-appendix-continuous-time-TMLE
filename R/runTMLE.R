### runTMLE.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  6 2020 (11:32) 
## Version: 
## Last-Updated: Nov  6 2020 (14:01) 
##           By: Thomas Alexander Gerds
##     Update #: 20
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
                    run.ctmle2 = TRUE,     # run continuous tmle 
                    misspecify.Q = FALSE,  # if TRUE, then misspecified model for outcome
                    only.A0 = FALSE,       # only baseline effect (not interesting here)
                    M = 5, # 500           # number of simulation repetitions
                    n = 1000,              # sample size
                    B = 5,                 # number of simulations to calculate true parameter value
                    N = 1000,              # sample size to calculate true parameter value
                    run.ctmle= FALSE,
                    seed,
                    progress.bar=4){

    #-------------------------------------------------------------------------------------------#
    ## true values (outputs to file)
    #-------------------------------------------------------------------------------------------#
    message("Computing true value of psi based on data generating mechanism:\n")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- progress.bar
        pb <- txtProgressBar(max = B,
                             style = progress.bar,
                             width=20)
    }
    true <- foreach(m=1:B, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
                    ) %dopar% {
                        if (progress.bar>0) { setTxtProgressBar(pb, m)}
                        x <- compute.true(n=N,
                                          compute.true.psi=TRUE,
                                          compute.true.eic=TRUE,
                                          only.A0=only.A0,K=K,misspecify.Q=misspecify.Q,seed=seed+m)
                    }
    #-------------------------------------------------------------------------------------------#
    ## repeat simulations (parallelize)
    #-------------------------------------------------------------------------------------------#
    message("\nEstimating psi with TMLE based on observed data:\n")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 4
        pb <- txtProgressBar(max = M,
                             style = progress.bar,
                             width=20)
    }
    out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
                   ) %dopar% {
                       if (progress.bar>0) { setTxtProgressBar(pb, m)}
                       repeat.fun(m, K=K, n=n,
                                  only.A0=only.A0, run.ltmle=run.ltmle, run.ctmle=run.ctmle,
                                  run.ctmle2=run.ctmle2,
                                  misspecify.Q=misspecify.Q,seed=seed+m)
                   }
    if (progress.bar>0) {    cat("\n")}
    if (no_cores>1) stopImplicitCluster()
    res <- list(true=true,est=out)
    class(res) <- "watmle"
    res
}
######################################################################
### runTMLE.R ends here
