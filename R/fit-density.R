### fit-density.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiese Rytgaard & Thomas Alexander Gerds
## Created: Oct 22 2020 (16:13) 
## Version: 
## Last-Updated: Oct 22 2020 (16:14) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
fit.density <- function(dt, Avar, covars, subset=NULL, test=FALSE) {

    ## dt.covars <- setDT(expand.grid(lapply(covars, function(var) {
    ##     dt[, unique(.SD), .SDcols=var][[1]]
    ## })))

    ## names(dt.covars) <- covars

    if (length(subset)>0) {
        dt.tmp <- dt[get(subset)==1, c(Avar, covars), with=FALSE]
    } else {
        dt.tmp <- dt[, c(Avar, covars), with=FALSE]
    }
    
    dt.tmp[, id:=1:nrow(dt.tmp)]

    avals <- sort(unique(dt.tmp[, Avar, with=FALSE][[1]]))
    size.A <- length(avals)-1
    
    for (a in avals) {
        dt.tmp[, (paste0("Avar", a)):=1*(a<=get(Avar))]
    }

    A.melt <- melt(dt.tmp, id=c("id", Avar, covars))[order(id)]

    A.melt[, a:=numextract(variable)]
    A.melt <- A.melt[, -"variable"]

    ## form <- formula(paste0("Y~a+", paste0(covars, collapse="+")))

    ## #dt.covars <- unique(dt.tmp[, covars, with=FALSE])
    
    ## hazard.A <- matrix(1, nrow=nrow(unique(dt.covars)), ncol=size.A+1)#rep(1, size.A+1)
    ## for (a1 in 0:(size.A-1)) {
    ##     A.melt[, Y:=(get(Avar)==a1)]
    ##     dt.covars[, a:=a1]
    ##     fit.A <- glm(form, data=A.melt[value==1 & a>=a1], family=binomial())
    ##     hazard.A[,a1+1] <- predict(fit.A, newdata=dt.covars, type="response")
    ## }

    if (TRUE) {
        A.melt[, atrisk:=1*(get(Avar)>=a)]
        A.melt[, Y:=1*(get(Avar)==a)]
        A.melt[, a:=factor(a)]
        
        form <- paste0("Y~", paste0(covars, "*a", collapse="+"))
        fit.A <- glm(formula(form), data=A.melt[atrisk==1], family=binomial())

        dt.covars <- setDT(expand.grid(lapply(c(covars, Avar), function(var) {
            dt[, unique(.SD), .SDcols=var][[1]]
        })))

        names(dt.covars) <- c(covars, Avar)
        dt.covars[, a:=factor(get(Avar))]

        #for (a1 in 0:(size.A-1)) {
        dt.covars[, hazard.A:=predict(fit.A, newdata=dt.covars, type="response")]
        #}

        dt.covars <- dt.covars[order(L0, get(Avar))]
        dt.covars[, haz.cumprod:=cumprod(1-hazard.A), by=covars]
        dt.covars[, est:=c(1,haz.cumprod[-.N])*hazard.A, by=covars]

        setnames(dt.covars, "est", paste0("fit.", Avar))
        return(dt.covars[, -c("a", "hazard.A", "haz.cumprod")])
        
    } else {

        out <- do.call("rbind", lapply(0:size.A, function(a) {
            if (a==0) est <- hazard.A[, 1] else est <- apply(hazard.A[, 1:a, drop=FALSE], 1, function(xx) {
                prod(1-xx)
            })*hazard.A[, a+1]
            data.table(dt.covars[, -"a", with=FALSE], Avar=a, fit=est)
        }))

        setnames(out, c("Avar", "fit"), c(Avar, paste0("fit.", Avar)))
        return(out)  
    }
}


######################################################################
### fit-density.R ends here
