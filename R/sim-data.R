sim.data <- function(n, # sample size
                     K=1, # maximal length of followup
                     seed, # seed
                     censoring=TRUE, # if FALSE data are uncensored
                     intervention.A=NULL, # used to calculate the true parameter value
                     intervention.A0=NULL, # used to calculate the true parameter value 
                     form.dN.L=NULL, # how time of covariate monitoring depends on history through last value (Markov)
                     form.dN.A=NULL, # how time of treatment monitoring depends on baseline and current history through last value (Markov)
                     form.A0=NULL, # baseline treatment
                     form.A=NULL, # how treatment decision depends on history through last value (Markov)
                     form.Y=NULL, # how outcome depends on history through last value (Markov)
                     form.C=NULL, # censoring
                     form.L=NULL # covariate
                     ) {

    if (length(form.dN.L)==0) form.dN.L <- function(L0, dN.L.prev, L.prev, A.prev) -0.2-0.05*K-0.025*(K>7)-0.25*dN.L.prev-0.15*L0-0.1*(A.prev==1)+0.3*L.prev
    if (length(form.dN.A)==0) form.dN.A <- function(L0, dN.A.prev, L.prev, A.prev) -0.75-0.05*K-0.42*dN.A.prev+0.15*L0+0.3*(A.prev==2)+0.4*(A.prev==1)-0.25*L.prev
    if (length(form.C)==0) form.C <- function(L0, L.prev, A.prev, A0) -3.95+(K>40)*5-0.4*K^{2/3}-0.24*(K>2 & K<=4)-0.4*(K>4 & K<=9)-(K>9)*0.4*K^{1/5}+0.2*(K>25)*K^{1/4}+0.1*L0+0.2*(A0==1)+0.9*(A0==2)+2.15*L.prev
    if (length(form.L)==0) form.L <- function(L0, L.prev, A.prev, A0) 0.5-0.4*A0+0.15*L0-0.25*(A.prev==1)+0.4*L.prev
    if (length(form.A0)==0) form.A0 <- function(L0) (-0.1+0.25*L0)
    if (length(form.A)==0) form.A <- function(L0, L.prev, A.prev, A0) (-1+(1-A0)*0.6+(1-A.prev)*0.4+L.prev*0.6-0.15*(K>15)*L.prev)
    if (length(form.Y)==0) {
        form.Y <- function(L0, L.prev, A.prev, A0, dN.A.prev){-1.1- 0.33*K/3*(K>2 & K<=4)-0.25*K^{2/3}-0.25*(K>4 & K<=9)- (K>25 & K<45)*0.3*K^{1/5}- (K>75)*0.31+(K>85)*0.2-(K>25 & K<75)*0.5*K^{1/5}+0.6*(K>25)*K^{1/4} - 0.25*A.prev + 0.4*L.prev - 0.25*A0 + 0.35*L.prev*A0 + (K>75)*0.1*A0+(K>85)*0.01*A0}
    }
    if (length(seed)>0) {
        set.seed(seed)
    }
    
    rexpit <- function(x) rbinom(n=n, size=1, prob=plogis(x))
    rmulti <- function(x) {
        if (!is.matrix(x)) x <- matrix(x, ncol=1)
        out <- rep(0, n)
        for (p in 1:ncol(x)) {
            out <- p*rexpit(x[,p])*(out==0) + out 
        }
        out[out==0] <- ncol(x)+1
        return(out-1)
    }

    L0 <- sample(1:6, n, replace=1000)/6

    A0 <- rmulti(form.A0(L0))
    
    if (length(intervention.A)>0 & length(intervention.A0)==0) {
        A0 <- rmulti(intervention.A(L0, 0, 0, A0))
    }
    if (length(intervention.A0)>0){
        A0 <- rmulti(intervention.A0(L0, A0))
    }

    A.prev <- A0
    dN.A.prev <- rep(0, n)
    dN.L.prev <- rep(0, n)
    L.prev <- rep(0, n)
    Y.prev <- rep(0, n)
    C.prev <- rep(0, n)


    dt <- data.table(id=1:n, L0=L0, A0=A0)
    
    for (k in 1:K) {

        #-- generate outcome:
        
        Y1 <- Y.prev + rexpit(form.Y(L0, L.prev, A.prev, A0, dN.A.prev))*(1-C.prev)*(1-Y.prev)

        #-- generate covariates: 

        dN.L1 <- rexpit(form.dN.L(L0, dN.L.prev, L.prev, A.prev))
        L1 <- dN.L1*rexpit(form.L(L0, L.prev, A.prev, A0))+(1-dN.L1)*L.prev

        #-- generate treatment: 

        dN.A1 <- rexpit(form.dN.A(L0, dN.A.prev, L.prev, A.prev))

        
        A1 <- rmulti(form.A(L0, L.prev, A.prev, A0))
        
        if (length(intervention.A)>0){
            A1 <- rmulti(intervention.A(L0, L.prev, A.prev, A1))
        }
        A1 <- dN.A1*A1 + (1-dN.A1)*A.prev
        
        #-- generate censoring: 
        
        if (length(intervention.A)==0 & length(intervention.A0)==0 & censoring) {
            C1 <- (1-Y1)*(C.prev+(1-C.prev)*rexpit(form.C(L0, L.prev, A.prev, A0)))
        } else {
            C1 <- rep(0, n)
        }

        dt.tmp <- data.table(Y1, dN.L1, L1, dN.A1, A1, C1)
        names(dt.tmp) <- paste0(c("Y", "dN.L", "L", "dN.A", "A", "C"), k)
        
        dN.A.prev <- dN.A1
        A.prev <- A1
        dN.L.prev <- dN.L1
        L.prev <- L1
        Y.prev <- Y1
        C.prev <- C1
        
        dt <- cbind(dt, dt.tmp)

    }

    dt[, (paste0("Y", K+1)):=Y.prev + rexpit(form.Y(L0, L.prev, A.prev, A0, dN.A.prev))*(1-C.prev)*(1-Y.prev)]
    dt[]
}
