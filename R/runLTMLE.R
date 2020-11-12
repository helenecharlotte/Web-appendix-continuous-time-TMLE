runLTMLE <- function(K=10,
                     n = 1000,
                     misspecify.init = FALSE,
                     seed,
                     M = 5,
                     no_cores=1,
                     progress.bar=3,
                     ...){
    if (M==1) progress.bar <- -1
    message("\nEstimating psi with LTMLE based on observed data:\n")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
        pb <- txtProgressBar(max = M,style = progress.bar,width=20)
    }
    out <- foreach(m=1:M, .errorhandling="pass") %dopar% {
        if (progress.bar>0) { setTxtProgressBar(pb, m)}
        # generate data
        dt <- sim.data(n,seed=seed+m,censoring=TRUE,K=K,...)
        estLTMLE(dt)
    }
    if (progress.bar>0)cat("\n")
    class(out) <- "watmle"
    out
}
    
estLTMLE <- function(data){
    df.wide <- as.data.frame(cbind(data[, "L0"],
                                   dN.A0=rep(1, nrow(data)),
                                   data[, -c("L0", "id")]))
    col.order <- names(df.wide)[-1]
    Anodes <- col.order[substr(col.order, 1, 1)=="A"]
    Ynodes <- col.order[substr(col.order, 1, 1)=="Y"]
    Cnodes <- col.order[substr(col.order, 1, 1)=="C"]
    Nnodes <- col.order[substr(col.order, 1, 4)=="dN.A"]
    Lnodes <- setdiff(col.order, c("L0", Anodes, Ynodes, Cnodes))

        for (ii in Cnodes) {
            df.wide[, names(df.wide)==ii] <- BinaryToCensoring(is.censored=df.wide[, names(df.wide)==ii])
        }

        if (TRUE) {
            abar0 <- (df.wide[, Nnodes] == 1) * 0 + (df.wide[, Nnodes] == 0) * df.wide[, Anodes]
            abar1 <- (df.wide[, Nnodes] == 1) * 1 + (df.wide[, Nnodes] == 0) * df.wide[, Anodes]
            #abar <- (df.wide[, Nnodes] == 1) * 0 + (df.wide[, Nnodes] == 0) * 0
            
            #--- specify intervention step 2: when dA(dt)=1 : 
            det.g.fun <- function(data, current.node, nodes) {
                if (substr(names(data)[current.node], 1, 1)=="C") {
                    if (FALSE) {
                        C.node.index <- current.node
                        observed.C <- data[, C.node.index]
                        is.deterministic <- observed.C=="censored" | is.na(observed.C)
                        prob1 <- observed.C[is.deterministic]#observed.C[is.deterministic]
                        return(list(is.deterministic = is.deterministic, prob1 = prob1))
                    }
                    #return(NULL)
                } else {                
                    N.nodes.index <- grep("dN.A", names(data))
                    N.node.index <- max(N.nodes.index[N.nodes.index < current.node]) 
                    A.node.index <- current.node
                    N <- data[, N.node.index]
                    observed.A <- data[, A.node.index]
                    is.deterministic <- N == 0 | is.na(N)
                    prob1 <- observed.A[is.deterministic]
                    return(list(is.deterministic = is.deterministic, prob1 = prob1))
                }
            }

            r0 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             abar=as.matrix(abar0),
                                             #abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

            r1 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             abar=as.matrix(abar1),
                                             #abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))
        } else {

            r0 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             #abar=as.matrix(abar),
                                             abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             #deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

            r1 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             #abar=as.matrix(abar),
                                             abar=rep(1, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             #deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

        }
        
    if (is(r0, "try-error")) {
        ltmle.list.0 <- "ERROR"
    } else {
        ltmle.list.0 <- list(est=r0$estimates["tmle"],
                             sd=summary(r0)$treatment$std.dev)
    }
    if (is(r1, "try-error")) {
        ltmle.list.1 <- "ERROR"
    } else {
        ltmle.list.1 <- list(est=r1$estimates["tmle"],
                             sd=summary(r1)$treatment$std.dev)
    }
    A0 <- unlist(ltmle.list.0)
    A1 <- unlist(ltmle.list.1)
    names(A0) <- c("ltmle.A0","sd.A0")
    names(A1) <- c("ltmle.A1","sd.A1")
    return(c(A0, A1))
}
