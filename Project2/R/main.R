rm(list=ls())
library(parallel)
mc.cores <- 4

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
#lmc <- lmc[lmc$Class=="Cep_F",]

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

## fit all the I curve.
fit_I <- list()
class_I <- list()
p_act_I <- list()
for(i in 1:nrow(lmc_I)){
    f <- paste("../lmc/I/",lmc_I[i,1],".dat",sep="")
    lc_I <- read.table(f)
    K <- 3
    class_I[[i]] <- paste(lmc_I[i,2])
    omegas <- get_freqs(1,100,.1/diff(range(lc_I[,1])))
    rss <- mclapply(omegas,compute_rss,c(0),K,lc_I,mc.cores=mc.cores)
    p_act_I[[i]] <- lmc_I[i,3]
    X <- construct_design(2*pi/p_act_I[[i]],K,lc_I[,1])
    beta <- compute_params(2*pi/p_act_I[[i]],K,lc_I[,2],lc_I[,3]^{-2},X)
    fit_I[[i]] <- get_sinusoidal_params(beta)
}


## fit all the V curve.
fit_V <- list()
class_V <- list()
p_act_V <- list()
for (i in 1:nrow(lmc_V)){
    f <- paste("../lmc/V/",lmc_V[i,1],".dat",sep="")
    lc_V <- read.table(f)
    K <- 3
    class_V[[i]] <- paste(lmc_V[i,2])
    omegas <- get_freqs(1,100,.1/diff(range(lc_V[,1])))
    rss <- mclapply(omegas,compute_rss,c(0),K,lc_V,mc.cores=mc.cores)
    p_act_V[[i]] <- lmc_V[i,3]
    X <- construct_design(2*pi/p_act_V[[i]],K,lc_V[,1])
    beta <- compute_params(2*pi/p_act_V[[i]],K,lc_V[,2],lc_V[,3]^{-2},X)
    fit_V[[i]] <- get_sinusoidal_params(beta)
}

feature_matrix_I <- cbind(data.frame(matrix(unlist(fit_I), nrow=length(fit_I), byrow=T)),lmc_I[,3])
colnames(feature_matrix_I) <- c("beta0","amp1","amp2","amp3","rho1","rho2","rho3","p")

feature_matrix_V <- cbind(data.frame(matrix(unlist(fit_V), nrow=length(fit_V), byrow=T)),lmc_V[,3])
colnames(feature_matrix_V) <- c("beta0","amp1","amp2","amp3","rho1","rho2","rho3","p")

### PCA

PC_I = prcomp(feature_matrix_I,center=TRUE,scale.=TRUE)
PC_V = prcomp(feature_matrix_V,center=TRUE,scale.=TRUE)

summary(PC_I)
summary(PC_V)

print(PC_I)
print(PC_V)

save(list=ls(),"run.RData")

#library(rgl)
#period_ranks <- rank(unlist(p_act_I))
### phase align light curves, normalized magnitudes, construct data frame
#lc_I_shift <- list()
#for(i in 1:length(p_act_I)){
#    f <- paste("../lmc/I/",lmc_I[i,1],".dat",sep="")
#    lc_I <- read.table(f)
#    t <- ((lc_I[,1] + (lmc_I[i,3]*fit_I[[i]]$rho[1])/(2*pi)) %% lmc_I[i,3])/lmc_I[i,3]
#    m <- lc_I[,2] - mean(lc_I[,2])
#    lc_I_shift[[i]] <- cbind(rep(period_ranks[i],length(t)),t,m)
#}
#dat <- do.call(rbind,lc_I_shift)
#plot3d(dat[,1],dat[,2],dat[,3],alpha=0.02,xlab="Period Rank",ylab="Phase",zlab="Normalized Magnitude")
#writeWebGL(filename="cepheids")
#plot(PC_I$x[,4],PC_I$x[,1])
#
#col = rep(c("red","blue"),each=20)
#plot(PC_I$x[,1], PC_I$x[,2], pch="", main = "Your Plot Title", xlab = "PC 1", ylab = "PC 2")
#text(PC_I$x[,1], PC_I$x[,2], labels=rownames(PC_I$x), col = col)
