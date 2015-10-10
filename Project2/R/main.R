rm(list=ls())
#library(foreach)
#library(doMC)
#registerDoMC(4)

construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
}

get_freqs <-function(period_min,period_max,freq_del = 0.1/4000){
    freq_max <- 1/period_min
    freq_min <- 1/period_max
    return(2 * pi * seq(freq_min, freq_max, freq_del))
}

compute_params <- function(w,K,mag,weights,X){
    B <- t(X) %*% (X * weights)
    d <- t(X) %*% (mag * weights)
    return(solve(B,d))
}   

compute_rss <- function(w,K,lc){
    X <- construct_design(w,K,lc[,1])
    beta <- compute_params(w,K,lc[,2],lc[,3]^2,X)
    r <- (lc[,2] - X%*%beta)
    return(sum(lc[,3] * (r^2)))
}
get_sinusoidal_params <- function(beta){
    beta0 <- beta[1]
    amp <- rep(0,(length(beta)-1)/2)
    rho <- rep(0,(length(beta)-1)/2)
    for(ii in 1:length(amp)){
        amp[ii] <- sqrt(beta[2*ii]^2 + beta[2*ii + 1]^2)
        rho[ii] <- atan2(beta[2*ii],beta[2*ii+1])
    }
    return(list(beta0=beta0,amp=amp,rho=rho))
}

lmc <- read.table("../lmc.txt", header=TRUE)
lmc <- lmc[lmc$Class=="Cep_F",]

## extract data for all the Cepheids
f_exists_I <- rep(TRUE,nrow(lmc))
f_exists_V <- rep(TRUE,nrow(lmc))

for(i in 1:nrow(lmc)){
    f_I <- paste("../lmc/I/",lmc[i,1],".dat",sep="")
    f_V <- paste("../lmc/V/",lmc[i,1],".dat",sep="")
    if(!file.exists(f_I)){
        f_exists_I[i] <- FALSE
    }
    if(!file.exists(f_V)){
        f_exists_V[i] <- FALSE
    }
} 
lmc_I <- lmc[f_exists_I,]
lmc_V <- lmc[f_exists_V,]

## fit all the Cepheids I curve.
fit_I <- list()
#lc_I <- list()
p_act_I <- list()
for(i in 1:nrow(lmc_I)){
    f <- paste("../lmc/I/",lmc_I[i,1],".dat",sep="")
    lc_I <- read.table(f)
    K <- 3
    omegas <- get_freqs(1,100,.1/diff(range(lc_I[,1])))
    rss <- vapply(omegas,compute_rss,c(0),K,lc_I)
    p_act_I[[i]] <- lmc_I[i,3]
    X <- construct_design(2*pi/p_act_I[[i]],K,lc_I[,1])
    beta <- compute_params(2*pi/p_act_I[[i]],K,lc_I[,2],lc_I[,3]^{-2},X)
    fit_I[[i]] <- get_sinusoidal_params(beta)
}


## fit all the Cepheids V curve.
fit_V <- list()
#lc_V <- list()
p_act_V <- list()
for (i in 1:nrow(lmc_V)) %dopar% {
    f <- paste("../lmc/V/",lmc_V[i,1],".dat",sep="")
    lc_V <- read.table(f)
    K <- 3
    omegas <- get_freqs(1,100,.1/diff(range(lc_V[,1])))
    rss <- vapply(omegas,compute_rss,c(0),K,lc_V)
    p_act_V[[i]] <- lmc_V[i,3]
    X <- construct_design(2*pi/p_act_V[[i]],K,lc_V[,1])
    beta <- compute_params(2*pi/p_act_V[[i]],K,lc_V[,2],lc_V[,3]^{-2},X)
    fit_V[[i]] <- get_sinusoidal_params(beta)
}
