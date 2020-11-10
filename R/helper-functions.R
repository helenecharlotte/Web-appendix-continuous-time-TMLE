print.watmle <- function(x,...){
    # estimated values
    x <- data.table(do.call("rbind",x))
    print(x[])
    x[]
}

#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#
summary.watmle <- function(object,true){
    x <- data.table(do.call("rbind",object))
    print(true)
    print(x)
    x[,{
        bias.A0 <- mean(psi.hat.A0)-true[["psi0.A0"]]
        bias.A1 <- mean(psi.hat.A1)-true[["psi0.A1"]]
        bias.diff <- (psi.hat.A1-psi.hat.A0)-(true[["psi0.A1"]]-true[["psi0.A0"]])
        cov.A0 <- 
    browser()
}


summary.watmle.old <- function(object,...){
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



