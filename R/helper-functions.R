print.watmle <- function(x,...){
    # estimated values
    X <- data.table(do.call("rbind",x))
    print(X[])
    X[]
}

#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

summary.watmle <- function(object,true){
    x <- data.table(do.call("rbind",object))
    corl <- length(grep("ltmle",names(x))>0)
    if (corl){
        setnames(x,gsub("ltmle","est",names(x)))
        setnames(x,gsub("sd","se",names(x)))
    }
    else{
        setnames(x,gsub("ctmle","est",names(x)))
    }
    out <- x[,{
        mean.A0 <- mean(est.A0)
        mean.A1 <- mean(est.A1)
        bias.A0 <- mean.A0-true[["psi0.A0"]]
        bias.A1 <- mean.A1-true[["psi0.A1"]]
        diff <- est.A1-est.A0
        true.diff <- true[["psi0.A1"]]-true[["psi0.A0"]]
        bias.diff <- mean(diff)-true.diff
        cov.A0 <- cov.fun(est=est.A0,se=se.A0,true=true[["psi0.A0"]])
        cov.A1 <- cov.fun(est=est.A1,se=se.A1,true=true[["psi0.A1"]])
        cov.diff <- cov.fun(est=diff,se=sqrt(se.A1^2+se.A1^2),true=true.diff)
        se.A0 <- mean(se.A0)
        se.A1 <- mean(se.A1)
        se.diff <- mean(sqrt(se.A0^2+se.A1^2))
        mse.A0 <- mse(est.A0)
        mse.A1 <- mse(est.A1)
        mse.diff <- mse(diff)
        out <- data.frame(matrix(c(true.A0=true[["psi0.A0"]],true.A1=true[["psi0.A1"]],true.diff=true.diff,
                                   mean.A0=mean.A0,mean.A1=mean.A1,mean.diff=mean.A1-mean.A0,
                                   bias.A0=bias.A0,bias.A1=bias.A1,bias.diff=bias.diff,
                                   se.A0=se.A0,se.A1=se.A1,se.diff=se.diff,
                                   cov.A0=cov.A0,cov.A1=cov.A1,cov.diff=cov.diff,
                                   mse.A0=mse.A0,mse.A1=mse.A1,mse.diff=mse.diff),ncol=3,byrow=TRUE))
        names(out) <- c("A0","A1","psi")
        out <- cbind(Result=c("true","mean","bias","se","coverage","MSE"),out)
        setnames(out,"Result",ifelse(corl,"LTMLE","cTMLE"))
    }]
    out
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

cov.fun <- function(est, se, true) {
    return(mean(est + qnorm(0.025)*se <= true &
                est + qnorm(0.975)*se >= true))
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



