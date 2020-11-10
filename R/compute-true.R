#-------------------------------------------------------------------------------------------#
## true values of parameters
#-------------------------------------------------------------------------------------------#

compute.true.psi0 <- function(n,
                              B=1,
                              K,
                              misspecify.init=FALSE,
                              seed,
                              verbose=TRUE,
                              no_cores=1,
                              progress.bar=TRUE) {
    stopifnot(K>0)
    stopifnot(is.numeric(seed))
    if (verbose==TRUE) message("Computing true value of psi based on data generating mechanism:\n")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {pb <- txtProgressBar(max = B,style = 3,width=20) }
    true <- foreach(m=1:B, .errorhandling="pass",.combine="rbind") %dopar% {
        if (progress.bar>0) { setTxtProgressBar(pb, m)}
        dt.A0 <- sim.data(n=n, seed=seed+m,
                          intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(1)),
                          K=K
                          )
        psi0.A0 <- mean(dt.A0[[paste0("Y", K+1)]])
        dt.A1 <- sim.data(n=n, seed=seed+B+m+1,
                          intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(0)),
                          K=K
                          )
        psi0.A1 <- mean(dt.A1[[paste0("Y", K+1)]])
        c(psi0.A0=psi0.A0,psi0.A1=psi0.A1)
    }
    if (progress.bar>0) { cat("\n")}
    if (B>1)
        psi0 <- apply(true,2,mean)
    else
        psi0 <- true
    if (verbose>1) {
        if (B<10)
            print(true)
        else
            print(apply(true,2,summary))
    }
    psi0
}
