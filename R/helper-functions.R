runTMLE <- function(no_cores,
                    K=10,#100              # end of follow-up (tau)
                    run.ltmle = FALSE,     # run ltmle, OR (!!)
                    run.ctmle2 = TRUE,     # run continuous tmle 
                    misspecify.Q = FALSE,  # if TRUE, then misspecified model for outcome
                    only.A0 = FALSE,       # only baseline effect (not interesting here)
                    M = 5, # 500           # number of simulation repetitions
                    n = 1000,              # sample size
                    run.ctmle= FALSE,
                    progress.bar=TRUE){

    #-------------------------------------------------------------------------------------------#
    ## true values (outputs to file)
    #-------------------------------------------------------------------------------------------#
    message("Computing true value of psi based on data generating mechanism")
    x <- compute.true(n=1000,
                      compute.true.psi=TRUE,
                      compute.true.eic=TRUE,
                      only.A0=only.A0)

    #-------------------------------------------------------------------------------------------#
    ## repeat simulations (parallelize)
    #-------------------------------------------------------------------------------------------#
    message("Estimating psi with TMLE based on observed data")
    if (no_cores>1) registerDoParallel(no_cores)
    if (progress.bar>0) {
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
        pb <- txtProgressBar(max = M,
                             style = progress.bar,
                             width=20)
    }
    out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
                   ) %dopar% {
                       setTxtProgressBar(pb, m)
                       repeat.fun(m, K=K, n=n,
                                  only.A0=only.A0, run.ltmle=run.ltmle, run.ctmle=run.ctmle,
                                  run.ctmle2=run.ctmle2,
                                  misspecify.Q=misspecify.Q)
                   }

    if (no_cores>1) stopImplicitCluster()
    res <- list(true=x,est=out)
    class(res) <- "watmle"
    res
}

print.watmle <- function(x,...){
    # true values
    psi0.A0 <- x$true$psi0.test.multi.M0
    psi0.A1 <- x$true$psi0.test.multi.M1
    # estimated values
    mat <- do.call("rbind",lapply(x$est,function(e){
        if (e[1]!="ERROR"){
            n1 <- length(e[[1]])
            n2 <- length(e[[2]])
            data.table(psiA0=e[[2]][[n2]][1],
                       psiA1=e[[1]][[n1]][1],
                       sdA0=e[[2]][[1]][3],
                       sdA1=e[[1]][[1]][3])
        }
    }))
    mat[,true.psiA0:=psi0.A0]
    mat[,true.psiA1:=psi0.A1]
    print(mat[])
    mat[]
}

#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

summary.watmle <- function(object,...){
    capture.output(X <- print(object))
    with(X,{
        cat("----------------------------","\n")
        cat("estimates from tmle:","\n") 
        cat("------------------","\n")
        cat(paste0("init (A=0): ", round(mean(psiA0.init), 4)),"\n")
        cat(paste0("tmle (A=0): ", round(mean(psiA0), 4)),"\n")
        cat(paste0("true (A=0): ", round(true.psiA0, 4)),"\n")
        cat("------------------","\n")
        cat(paste0("init (A=1): ", round(mean(psiA1.init), 4)),"\n")
        cat(paste0("tmle (A=1): ", round(mean(psiA1), 4)),"\n")
        cat(paste0("true (A=1): ", round(true.psiA1, 4)),"\n")
        cat("-----------------------------------","\n")
        cat("coverage of tmle estimator:","\n")
        cat(paste0("coverage (A=0): ", round(cov.fun(psiA0, sdA0, true.psiA0), 4)),"\n")
        cat(paste0("coverage (A=1): ", round(cov.fun(psiA1, sdA1, true.psiA1), 4)),"\n")
        cat(paste0("coverage (diff): ", round(cov.fun(psiA1-psiA0, sqrt(sdA1^2+sdA0^2), true.psiA1-true.psiA0), 4)),"\n")
        cat("----------------------------","\n")
        cat("se estimates (A=0):","\n")
        cat(paste0("mean sigma : ", round(mean(sdA0), 4)),"\n")
        cat(paste0("mse        : ", round(mse(psiA0), 4)),"\n")
        cat("----------------------------","\n")
        cat("se estimates (A=1):","\n")
        cat(paste0("mean sigma : ", round(mean(sdA1), 4)),"\n")
        cat(paste0("mse        : ", round(mse(psiA1), 4)),"\n")
        cat("-------------------------------","\n")
        cat("efficiency (A=0):","\n")
        cat(paste0("MSE / sd       : ", round(mse(psiA0)/mean(sdA0), 4)),"\n")
        cat("-------------------------------","\n")
        cat("efficiency (A=1):","\n")
        cat(paste0("MSE / sd       : ", round(mse(psiA1)/mean(sdA1), 4)),"\n")
    })
}

plot.watmle <- function(x,...){
    X <- print(x)
    with(X,{
        hist(sdA1)
        hist(sdA0)
        par(mfrow=c(2,2))
        hist(psiA0.init)
        abline(v=true.psiA0, col="red")
        hist(psiA1.init)
        abline(v=true.psiA1, col="red")        
        hist(psiA0)
        abline(v=true.psiA0, col="red")
        hist(psiA1)
        abline(v=true.psiA1, col="red")
    })
}

cov.fun <- function(psi.hat, sd, psi0) {
    return(mean(psi.hat - 1.96*sd <= psi0 &
                psi.hat + 1.96*sd >= psi0))
}

numextract <- function(string){ 
    as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

logit <- function(p) log(p/(1-p)) #qlogis :) 
mse <- function(x, x0=NULL) {
    if (length(x0)==0) {
        return(sqrt(mean((x-mean(x))^2)))
    } else {
        return(sqrt(mean((x-x0)^2)))
    }
}



