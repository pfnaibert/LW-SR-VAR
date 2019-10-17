##########################
# Paulo Ferreira Naibert
# ORIGINAL: 13 of October of 2019
# LAST MODIFIED: 17 of October of 2019
##########################

##########################
diff.test <- function(ret, est, rav=1, kernel="parzen", pw=0) 
{

if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

if(!(est %in% c("var", "sr", "mu", "ceq"))){stop("\n est must be either 'var', 'sr', 'mu', or 'ceq' \n\n")}

if(!(kernel %in% c("parzen", "sample"))){stop("\n kernel must be either 'sample' or 'parzen' \n\n")}

obj <- est.diff(ret, est, rav); D.hat <- obj$Diff;

se <- se.fn(ret, est, rav, kernel, pw);
test  <- -abs(D.hat)/se; PV <- 2*pnorm(test);

out <- list("Diff"=D.hat, "se"=se, "test"=test, "p.value"=PV)
return(out)
}

##########################
se.fn <- function(ret, est, rav=1, kernel="parzen", pw=0) 
{
if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

if(!(est %in% c("var", "sr", "mu", "ceq"))){stop("\n est must be either 'var', 'sr', 'mu', or 'ceq' \n\n")} else

if(!(kernel %in% c("parzen", "sample"))){stop("\n kernel must be either 'sample' or 'parzen' \n\n")}

obj <- util.fn(ret, est, rav)
T   <-  NROW(ret); gradient <- obj$grad; V.hat <- obj$V.hat

if(kernel=="sample"){Psi.hat <- cov(V.hat)} else

if(kernel=="parzen")
{
if(pw==0){Psi.hat <- Psi.hat.fn(V.hat)}

if(pw==1)
{
tmp     <- prewhite.fn(V.hat); V.star  <- tmp$V.star; D <- tmp$D
Psi.hat <- Psi.hat.fn(V.star); Psi.hat <- D%*%tcrossprod(Psi.hat, D)
}
}

return(as.numeric(sqrt((1/T)*crossprod(gradient, Psi.hat%*%gradient))) )
}

##########################
Psi.hat.fn <- function(V.hat) 
{
T <- NROW(V.hat)

alpha.hat <- alpha.hat.fn(V.hat)
S.star    <- 2.6614*(alpha.hat*T)^0.2; S.star <- min(S.star, T-1)
Psi.hat   <- Gamma.hat.fn(V.hat, 0)

j <- 1
while (j < S.star)
{
Gamma.hat <- Gamma.hat.fn(V.hat, j)
Psi.hat   <- Psi.hat + kernel.Parzen(j/S.star)*(Gamma.hat + t(Gamma.hat))
j         <- j + 1
}

return((T/(T-4))*Psi.hat)
}

##########################
Gamma.hat.fn <- function(V.hat, j) 
{
T <- NROW(V.hat); p <- NCOL(V.hat); Gamma.hat <- matrix(0, p, p)
if(j >= T){ stop("j must be smaller than the row dimension!") }
for(i in ((j + 1):T)){Gamma.hat <- Gamma.hat + tcrossprod(V.hat[i, ], V.hat[i - j, ]) } 

return(Gamma.hat/T)
}

##########################
alpha.hat.fn <- function(V.hat) 
{

T <- NROW(V.hat); p <- NCOL(V.hat); num <- 0; den <- 0

for(i in (1:p))
{
fit <- ar(V.hat[, i], 0, 1, method = "ols")
rho.hat <- as.numeric(fit[2])
sig.hat <- sqrt(as.numeric(fit[3]))
num <- num + 4*rho.hat^2*sig.hat^4/(1 - rho.hat)^8
den <- den + sig.hat^4/(1 - rho.hat)^4
}

return(num/den)
}

##########################
kernel.Parzen <- function(x) 
{
if(abs(x) <= 0.5){result <- 1 - 6*x^2 + 6*abs(x)^3}
else if(abs(x) <= 1){result <- 2*(1 - abs(x))^3}
else result <- 0
return(result)
}

##########################
# Handles V.hat with N columns
# PW handles only VAR models
prewhite.fn <- function(V.hat)
{

T <- NROW(V.hat); N <- NCOL(V.hat)
seq1 <- seq(1, T-1); seq2 <- seq(2, T)

A.ls <- matrix(0, N, N); V.star <- matrix(0, T-1, N)        
vars <- matrix(0, T-1, N)

for(i in (1:N)){vars[,i] <- V.hat[seq1, i]}

# VAR FIT
for(j in (1:N))
{
fit      <- lm(V.hat[seq2, j] ~ -1 + vars)
A.ls[j,] <- as.numeric(fit$coef); V.star[, j] <- as.numeric(fit$resid)
}

svd.A <- svd(A.ls); d <- svd.A$d; d.adj <- d

# VAR COEFS ADJ
for(i in (1:N))
{
if(d[i] >  0.97){d.adj[i] <- 0.97} else
if(d[i] < -0.97){d.adj[i] <- -0.97}
}

A.hat   <- svd.A$u%*%tcrossprod(diag(d.adj), svd.A$v)
D       <- solve(diag(N) - A.hat)
reg.mat <- t(vars)

for(j in (1:N)){V.star[,j] <- V.hat[seq2,j] - A.hat[j,]%*%reg.mat}

out <- list("V.star"=V.star, "D"=D)
return(out)
}

#####################################
boot.test <- function(ret, est, rav=1, b=5, M=499, D.null=0) 
{
if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

if(!(est %in% c("var", "sr", "mu", "ceq"))){stop("\n est must be either 'var', 'sr', 'mu', or 'ceq' \n\n")} else

# 
T <- NROW(ret); l <- floor(T/b);
D.hat <- est.diff(ret, est, rav)$Diff;
d <- abs(D.hat-D.null)/se.fn(ret, est, rav, kernel="parzen", pw=1);
p.value <- 1

#
cat("\n ========================== ")
for(m in (1:M))
{

#
if(m%%500==0){cat("\n *** Running Bootstrap with block length", b, "Iteration", m, "out of", M, "***")}

#
ret.star   <- ret[cbb.seq(T,b),];
D.hat.star <- est.diff(ret.star, est, rav)$Diff; 

obj <- util.fn(ret.star, est, rav); gradient <- obj$grad; y.star <- obj$V.hat
N   <- NCOL(y.star)

Psi.hat.star <- matrix(0, N, N)
for(j in (1:l))
{
zeta.star <- sqrt(b)*colMeans(y.star[(1 + (j-1)*b):(j*b),])
Psi.hat.star <- Psi.hat.star + tcrossprod(zeta.star)
}
Psi.hat.star <- (T/(T-4))*Psi.hat.star/l

se.star <- as.numeric(sqrt((1/T)*crossprod(gradient, Psi.hat.star%*%gradient)))
d.star  <- abs(D.hat.star - D.hat)/se.star

if(d.star >= d){p.value <- p.value + 1}
}
cat("\n ========================== \n")

p.value <- p.value/(M + 1)

out <- list("Diff"=D.hat, "p.value"=p.value)
return(out)
}

#####################################
cbb.seq <- function(T, b)
{
# Circluar bootstrap
# Boot data is made of l blocks of size b
# circular bootstrap has UNIT MASS (only b=5 allowed)

T1   <- seq(1, T); l <- floor(T/b);
ids0 <- c(T1, 1:b)    # wrap the data around a circle # original indexes
ids1 <- rep(0, T)  # bootstrap indexes 
sps  <- sample(T1, l,replace=TRUE) # index where the block starts

for(j in seq(1, l))
{
ids1[seq((1 + (j-1)*b), (j*b))] <- ids0[seq(sps[j], (sps[j] + b - 1))]
}

return(ids1)
}

#####################################
sb.seq <- function(T, b.av) 
{
# Stationary bootstrap
# Boot data is made of l block with average size b

id.seq <- rep(1:T, 2) # wrapping data
seq1   <- rep(NA, 2*T)  # pre-allocate seq vector
start  <- sample(1:T, T, replace = T)
b      <- 1 + rgeom(T, 1/b.av)            # took the RNG out of the loop

i=0
while (i < T)
{
seq1[seq(1+i, (i+b[1+i]) )] <- id.seq[start[1+i]:(start[1+i] + b[1+i] - 1)]
i <- i + b[1+i]
}

out <- seq1[1:T]

return(out)
}

##########################
est.diff <- function(ret, est, rav=1)
{
if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

if(!(est %in% c("var", "sr", "mu", "ceq"))){stop("\n est must be either 'var', 'sr', 'mu', or 'ceq' \n\n")} else

if(est=="mu"){est1  <-  mean(ret[,1]); est2 <-  mean(ret[,2])} else
if(est=="var"){est1 <-  log(var(ret[,1])); est2 <-  log(var(ret[,2]))} else
if(est=="ceq"){est1 <- mean(ret[,1]) - (rav/2)*var(ret[,1]); est2 <- mean(ret[,2]) - (rav/2)*var(ret[,2])} else
if(est=="sr"){est1  <- mean(ret[,1])/sd(ret[,1]); est2 <- mean(ret[,2])/sd(ret[,2])} 

Diff <- est1 - est2

out <- list("Diff"=Diff, "ests"=c("est1"=est1, "est1"=est2))
return(out)
}

#####################################
# utility function to calculate the V.hat AND the gradient
util.fn <- function(ret, est, rav=1)
{

if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

if(!(est %in% c("var", "sr", "mu", "ceq"))){stop("\n est must be either 'var', 'sr', 'mu', or 'ceq' \n\n")}

dat   <- cbind(ret, ret^2); mu <- apply(dat, 2, mean)
V.hat <- cbind(dat[,1]-mu[1], dat[,2]-mu[2], dat[,3]-mu[3], dat[,4]-mu[4])
gradient <- rep(NA, 4)

if(est=="mu")
{
V.hat    <- V.hat[,c(1,2)]
gradient <- c(1, -1)
}

if(est=="sr")
{
# sr grad
gradient[1] <-  mu[3]/(mu[3]-mu[1]^2)^1.5
gradient[2] <- -mu[4]/(mu[4]-mu[2]^2)^1.5
gradient[3] <- -(0.5)*mu[1]/(mu[3]-mu[1]^2)^1.5
gradient[4] <-  (0.5)*mu[2]/(mu[4]-mu[2]^2)^1.5
}
   
if(est=="var")
{
# var grad
gradient[1] <- -2*mu[1]/(mu[3] - mu[1]^2)
gradient[2] <-  2*mu[2]/(mu[4] - mu[2]^2)
gradient[3] <-  1/(mu[3] - mu[1]^2)
gradient[4] <- -1/(mu[4] - mu[2]^2)
}

if(est=="ceq")
{
# ceq grad
gradient[1] <-  1 + rav*mu[1]
gradient[2] <- -1 - rav*mu[2]
gradient[3] <- -rav/2
gradient[4] <-  rav/2
}

out <- list("grad"=gradient, "V.hat"=V.hat)
return(out)
}

#####################################
