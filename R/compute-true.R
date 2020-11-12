#-------------------------------------------------------------------------------------------#
## true values of parameters
#-------------------------------------------------------------------------------------------#

compute.true <- function(K, # time grid
                         n, # sample size
                         B=1, # number of repetitions 
                         seed,
                         no_cores=1,
                         progress.bar=TRUE,...) {
    stopifnot(K>0)
    stopifnot(is.numeric(seed))
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {pb <- txtProgressBar(max = B,style = 3,width=20) }
    true <- foreach(m=1:B, .errorhandling="pass",.combine="rbind") %dopar% {
        if (progress.bar>0) { setTxtProgressBar(pb, m)}
        dt.A0 <- sim.data(n=n, seed=seed+m,
                          intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(1)),
                          K=K,...)
        psi0.A0 <- mean(dt.A0[[paste0("Y", K+1)]])
        dt.A1 <- sim.data(n=n, seed=seed+B+m+1,
                          intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(0)),
                          K=K,...)
        Y.A1 <- dt.A1[[paste0("Y", K+1)]]
        Y.A0 <- dt.A0[[paste0("Y", K+1)]]
        cbind(Y.A0,Y.A1)
    }
    psi0 <- apply(true,2,mean)
    names(psi0) <- c("psi0.A0","psi0.A1")
    if (progress.bar>0) { cat("\n")}
    psi0
}
