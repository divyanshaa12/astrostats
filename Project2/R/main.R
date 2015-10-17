rm(list=ls())
#library(rainbow)
#library(parallel)
#mc.cores <- 1

construct_design <- function(w,K,t){
    predesign <- w*outer(t,1:K)
    return(cbind(1,cos(predesign),sin(predesign)))
}

get_freqs <-function(period.min,period.max,freq.del = 0.1/4000){
    freq.max <- 1/period.min
    freq.min <- 1/period.max
    return(2 * pi * seq(freq.min, freq.max, freq.del))
}

compute_params <- function(w,K,mag,weights,X){
    B <- t(X) %*% (X * weights)
    d <- t(X) %*% (mag * weights)
    return(solve(B,d))
}   

compute_rss <- function(w,K,lc){
    X <- construct.design(w,K,lc[,1])
    beta <- compute.params(w,K,lc[,2],lc[,3]^2,X)
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

lmc.all <- read.table("../lmc.txt", header=TRUE)

# limiting number of light curves for each class to a max of 10. Just testing
classes <- paste(unique(lmc.all[,2]))
rm(out,tmp)
tmp <- lmc.all[lmc.all$Class == classes[1],]
out <- tmp[1:min(10,nrow(tmp)),]
for (i in 2:length(classes)){
    tmp <- lmc.all[lmc.all$Class == classes[i],]
    tmp <- tmp[1:min(10,nrow(tmp)),]
    out <- rbind(out,tmp)
}

lmc <- out

f.exists.I <- rep(TRUE,nrow(lmc))
f.exists.V <- rep(TRUE,nrow(lmc))

for(i in 1:nrow(lmc)){
    f.I <- paste("../lmc/I/",lmc[i,1],".dat",sep="")
    f.V <- paste("../lmc/V/",lmc[i,1],".dat",sep="")
    if(!file.exists(f.I)){
        f.exists.I[i] <- FALSE
    }
    if(!file.exists(f.V)){
        f.exists.V[i] <- FALSE
    }
} 
lmc.I <- lmc[f.exists.I,]
lmc.V <- lmc[f.exists.V,]
lmc.I <- lmc.I[lmc.I$Period >0,]
lmc.V <- lmc.V[lmc.V$Period >0,]



M=49
## fit all the I curve.
fit.std.I <-  matrix(0,M+1,nrow(lmc.I))
p.act.I <- list()
for(i in 1:nrow(lmc.I)){
  f <- paste("../lmc/I/",lmc.I[i,1],".dat",sep="")
  lc.I <- read.table(f)
  K <- 8
  omegas <- get_freqs(1,100,.1/diff(range(lc.I[,1])))
  p.act.I[[i]] <- lmc.I[i,3]
  X <- construct_design(2*pi/p.act.I[[i]],K,lc.I[,1])
  beta <- compute_params(2*pi/p.act.I[[i]],K,lc.I[,2],lc.I[,3]^{-2},X)
  
  grid = seq(min(lc.I[,1]%%p.act.I[[i]]),max(lc.I[,1]%%p.act.I[[i]]),(max(lc.I[,1]%%p.act.I[[i]])-min(lc.I[,1]%%p.act.I[[i]]))/M)
  X.grid = construct_design(2*pi/p.act.I[[i]],K,grid)
  y.grid = X.grid%*%beta
  fit.std.I[,i] = y.grid
  if(i%%500==0){print(i)}
}


## fit all the V curve.
 fit.std.V <-  matrix(0,M+1,nrow(lmc.V))
p.act.V <- list()
for(i in 1:nrow(lmc.V)){
  f <- paste("../lmc/I/",lmc.V[i,1],".dat",sep="")
  lc.V <- read.table(f)
  K <- 8
  omegas <- get_freqs(1,100,.1/diff(range(lc.V[,1])))
  p.act.V[[i]] <- lmc.V[i,3]
  X <- construct_design(2*pi/p.act.V[[i]],K,lc.V[,1])
  beta <- compute_params(2*pi/p.act.V[[i]],K,lc.V[,2],lc.V[,3]^{-2},X)
  
  grid = seq(min(lc.V[,1]%%p.act.V[[i]]),max(lc.V[,1]%%p.act.V[[i]]),(max(lc.V[,1]%%p.act.V[[i]])-min(lc.V[,1]%%p.act.V[[i]]))/M)
  X.grid = construct_design(2*pi/p.act.V[[i]],K,grid)
  y.grid = X.grid%*%beta
  fit.std.V[,i] = y.grid
  if(i%%500==0){print(i)}
}

Fa.I = cov(t(fit.std.I))
Fa.I.e = eigen(Fa.I)
Fa.I.ve = Fa.I.e$vectors
Fa.I.va = Fa.I.e$values

Fa.V = cov(t(fit.std.V))
Fa.V.e = eigen(Fa.V)
Fa.V.ve = Fa.V.e$vectors
Fa.V.va = Fa.V.e$values

H=3
n = ncol(fit.std.I) 
coef = matrix(0,n,H)
for(i in 1:n){
  fit1 = lm(fit.std.I[,i]~0+Fa.I.ve[,1:H])
  coef[i,] = fit1$coefficients
  if(i%%500==0){print(i)}
}

ind.I <- list()
for (i in 1:length(classes)){
    ind.I[[i]] <- which(lmc.I[,2]==classes[i])
}
ind.I <- unlist(ind.I)

#
####### Scatterplot
#colnames(coef) = c("PC1","PC2","PC3")
#dat$class <- as.factor(as.character(lmc_I[ind0,2]))
#dat$feature = cbind(coef[ind0,],lmc_I$Period[ind0])
#pairs(coef[ind0,1:H],col = c("orange","blue","black","red","green")[dat$class],pch=(1:6)[dat$class])
#
####### 3D plot
#library(rgl)
#plot3d(coef[ind1,1], coef[ind1,2], coef[ind1,3], col="green", size=2)
#plot3d(coef[ind2,1], coef[ind2,2], coef[ind2,3], col="blue", size=2,add=TRUE)
#plot3d(coef[ind3,1], coef[ind3,2], coef[ind3,3], col="red", size=2,add=TRUE)
#plot3d(coef[ind4,1], coef[ind4,2], coef[ind4,3], col="black", size=2,add=TRUE)
#plot3d(coef[ind5,1], coef[ind5,2], coef[ind5,3], col="orange", size=2,add=TRUE)
#
#plot3d(coef[ind1,1], coef[ind1,2], log(lmc_I$Period[ind1]), col="green", size=2)
#plot3d(coef[ind2,1], coef[ind2,2], log(lmc_I$Period[ind2]), col="blue", size=2,add=TRUE)
#plot3d(coef[ind3,1], coef[ind3,2], log(lmc_I$Period[ind3]), col="red", size=2,add=TRUE)
#plot3d(coef[ind4,1], coef[ind4,2], log(lmc_I$Period[ind4]), col="black", size=2,add=TRUE)
#plot3d(coef[ind5,1], coef[ind5,2], log(lmc_I$Period[ind5]), col="orange", size=2,add=TRUE)

###### Fit Random Forest

#n = length(ind)
#test.size = 10
#samp = sample(1:n,test.size)
#ind.te = ind[[1]][samp]
#ind.tr = ind[[1]][-samp]

dat.I = list()
dat.I$class <- as.factor(as.character(lmc.I[ind.I,2]))
dat.I$feature = cbind(coef[ind.I,],lmc.I$Period[ind.I])

#dat.te = list()
#dat.te$class <- as.factor(as.character(lmc.I[ind.te,2]))
#dat.te$feature = cbind(coef[ind.te,],lmc.I$Period[ind.te])

dat.I = data.frame(dat.I$feature,dat.I$class)
#dat.te = data.frame(dat.te$feature,dat.te$class)
colnames(dat.I) = c("PC1","PC2","PC3","Period","class")
#colnames(dat.te) = c("PC1","PC2","PC3","Period","class")

#colnames(dat.tr) = c("PC1","PC2","PC3","Period","class")
#colnames(dat.te) = c("PC1","PC2","PC3","Period","class")

rf.fit.I = rpart(class~.,data = dat.I)
#library(tables)
#tblr <- tabular((Truth=dat.te$class) ~ (Predicted = predict(fit,dat.te,type="class")))
#print(tblr)
#library(rpart.plot)
#prp(fit,extra=2,compress=FALSE,varlen=0)


n = ncol(fit.std.V) 
coef = matrix(0,n,H)
for(i in 1:n){
  fit1 = lm(fit.std.V[,i]~0+Fa.V.ve[,1:H])
  coef[i,] = fit1$coefficients
  if(i%%500==0){print(i)}
}

ind.V <- list()
for (i in 1:length(classes)){
    ind.V[[i]] <- which(lmc.V[,2]==classes[i])
}
ind.V <- unlist(ind.V)

dat.V = list()
dat.V$class <- as.factor(as.character(lmc.V[ind.V,2]))
dat.V$feature = cbind(coef[ind.V,],lmc.V$Period[ind.V])
dat.V = data.frame(dat.V$feature,dat.V$class)
colnames(dat.V) = c("PC1","PC2","PC3","Period","class")
rf.fit.V = rpart(class~.,data = dat.V)
