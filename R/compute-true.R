#-------------------------------------------------------------------------------------------#
## true values of parameters
#-------------------------------------------------------------------------------------------#

compute.true <- function(n=1e6,compute.true.psi=FALSE,compute.true.eic=FALSE,only.A0=FALSE,K,misspecify.Q,seed) {
    out <- NULL                         
    if (compute.true.psi) {
        
        psi0.test.multi.M0 <- sim.data(n=n, seed=seed,
                                       only.A0=only.A0,
                                       intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(1)),
                                       K=K
                                       )

        psi0.test.multi.M1 <- sim.data(n=n, seed=seed+1,
                                       only.A0=only.A0,
                                       intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(0)),
                                       K=K
                                       )
        out <- c(out,list(psi0.test.multi.M0=psi0.test.multi.M0,psi0.test.multi.M1=psi0.test.multi.M1))
    }

    if (compute.true.eic) {

        N <- pmin(n,1e5)
        dt <- sim.data(N, seed=seed, censoring=TRUE,
                       only.A0=only.A0,
                       browse=FALSE,
                       K=K)

        true.eic.0 <-  suppressMessages(est.fun(copy(dt), censoring=TRUE,
                                                targeting=1,
                                                only.A0=only.A0,
                                                smooth.initial=TRUE,
                                                browse9=FALSE, maxIter=1, 
                                                compute.true.eic=TRUE,
                                                intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                                intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                                                browse0=FALSE, misspecify.Q=misspecify.Q))
        
        true.eic.1 <- suppressMessages(est.fun(copy(dt), censoring=TRUE,
                                               targeting=1,
                                               only.A0=only.A0,
                                               smooth.initial=TRUE,
                                               browse9=FALSE, maxIter=1, 
                                               compute.true.eic=TRUE,
                                               intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                               intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==1)),
                                               browse0=FALSE, misspecify.Q=misspecify.Q))
        true.eic0 <- true.eic.0[[1]][3] * sqrt(N) / sqrt(1000)
        true.eic1 <- true.eic.1[[1]][3] * sqrt(N) / sqrt(1000)
        out <- c(out,list(true.eic0=true.eic0,true.eic1=true.eic1))
    }
} 
