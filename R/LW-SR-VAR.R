##########################
# Paulo Ferreira Naibert
# ORIGINAL: 13 of October of 2019
##########################

##########################
lv.diff <- function(ret)
{
lv1  <-  log(var(ret[,1])); lv2 <-  log(var(ret[,2])); Diff <- lv1 - lv2
out <- list("Diff"=Diff, "Vars"=c("Var1"=lv1, "Var2"=lv2))
return(out)
}

#####################################
sr.diff <- function(ret)
{
SR1 <- mean(ret[,1])/sd(ret[,1]); SR2 <- mean(ret[,2])/sd(ret[,2]); Diff <- SR1 - SR2;
out <- list("Diff"=Diff, "SRs"=c("SR1"=SR1, "SR2"=SR2))
return(out)
}

##########################
var.hac.test <- function(ret, pw=0) 
{
obj <- lv.diff(ret); D.hat <- obj$Diff

if(!(pw %in% c(0, 1))){stop("\n PW must be either '0' or '1' \n\n")} else
if(pw==0){se <- se.Parzen(ret, "var")} else
if(pw==1){se <- se.Parzen.pw(ret, "var")}

PV  <- 2*pnorm(-abs(D.hat)/se);

out <- list("Diff"=obj$Diff, "se"=se, "p.value"=PV)
return(out)
}

##########################
sr.hac.test <- function(ret, pw=0) 
{
obj <- sr.diff(ret); D.hat <- obj$Diff

if(!(pw %in% c(0, 1))){stop("\n PW must be either '0' or '1' \n\n")} else
if(pw==0){se <- se.Parzen(ret, "sr")}
if(pw==1){se <- se.Parzen.pw(ret, "sr")}
PV  <- 2*pnorm(-abs(D.hat)/se);

out <- list("Diff"=D.hat, "se"=se, "p.value"=PV)
return(out)
}

#####################################
sr.boot.1 <- function (ret, b, M, Delta.null = 0) 
{
T = NROW(ret); l = floor(T/b)

Delta.hat = sharpe.ratio.diff(ret)
d <- abs(Delta.hat-Delta.null)/se.Parzen.pw(ret, "sr"); # changed

p.value = 1

for (m in (1:M)) {
ret.star = ret[cbb.seq(T, b), ]; Delta.hat.star = sharpe.ratio.diff(ret.star)

obj <- sr.util(ret.star); gradient <- obj$grad; # y.star <- obj$V.hat

ret1.star = ret.star[, 1]
ret2.star = ret.star[, 2]
mu1.hat.star = mean(ret1.star)
mu2.hat.star = mean(ret2.star)
gamma1.hat.star = mean(ret1.star^2)
gamma2.hat.star = mean(ret2.star^2)

y.star = data.frame(ret1.star - mu1.hat.star, ret2.star - mu2.hat.star, ret1.star^2 - gamma1.hat.star, ret2.star^2 - gamma2.hat.star)

Psi.hat.star = matrix(0, 4, 4)
for (j in (1:l)) {
    zeta.star = sqrt(b)*colMeans(y.star[((j - 1) * b + 1):(j * 
        b), ])
    Psi.hat.star = Psi.hat.star + zeta.star %*% t(zeta.star)
}
Psi.hat.star = Psi.hat.star/l

Psi.hat.star = (T / (T - 4)) * Psi.hat.star
se.star = as.numeric(sqrt(t(gradient) %*% Psi.hat.star %*% gradient/T))
d.star = abs(Delta.hat.star - Delta.hat)/se.star

if(d.star >= d){p.value = p.value + 1}
}

p.value = p.value/(M + 1)
out <- list("Difference" = Delta.hat, "p.Value"=p.value )
return(out)
}

#####################################
sr.boot.test <- function(ret, b=5, M=499, D.null=0) 
{
if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

T <- NROW(ret); l <- floor(T/b);
D.hat <- sr.diff(ret)$Diff; d <- abs(D.hat-D.null)/se.Parzen.pw(ret, "sr");
p.value <- 1

cat("\n ========================== ")
for(m in (1:M))
{

if(m%%500==0){cat("\n *** Running Bootstrap with block length", b, "Iteration", m, "out of", M, "***")}

ret.star <- ret[cbb.seq(T,b),]; D.hat.star <- sr.diff(ret.star)$Diff; 
obj <- sr.util(ret.star); gradient <- obj$grad; y.star <- obj$V.hat

Psi.hat.star <- matrix(0, 4, 4)
for(j in (1:l))
{
zeta.star <- sqrt(b)*colMeans(y.star[(1 + (j-1)*b):(j*b),])
Psi.hat.star <- Psi.hat.star + tcrossprod(zeta.star)
}
Psi.hat.star <- (T/(T-4))*Psi.hat.star/l

se.star <- as.numeric(sqrt(crossprod(gradient, Psi.hat.star%*%gradient)/T))
d.star  <- abs(D.hat.star - D.hat)/se.star

if(d.star >= d){p.value <- p.value + 1}
}
cat("\n ========================== \n")

p.value <- p.value/(M + 1)

out <- list("Diff"=D.hat, "p.value"=p.value)
return(out)
}

##########################
var.boot.test <- function(ret, b=5, M=499, D.null = 0) 
{
if(NCOL(ret)!=2){stop("\n Number of series differ from 2 \n")}

T <- NROW(ret); l <- floor(T/b); 
D.hat <- lv.diff(ret)$Diff; d <- abs(D.hat - D.null)/se.Parzen.pw(ret, "var");
p.value <- 1

cat("\n ========================== ")
for(m in (1:M))
{

if(m%%500==0){cat("\n *** Running Bootstrap with block length", b, "Iteration", m, "out of", M, "***")}
ret.star <- ret[cbb.seq(T,b),]; D.hat.star <- lv.diff(ret.star)$Diff;
obj <- var.util(ret); gradient <- obj$grad; y.star <- obj$V.hat

Psi.hat.star <- matrix(0, 4, 4)
for(j in (1:l))
{
zeta.star <- sqrt(b)*colMeans(y.star[(1 + (j-1)*b):(j*b), ])
Psi.hat.star <- Psi.hat.star + tcrossprod(zeta.star)
}
Psi.hat.star <- (T/(T-4))*Psi.hat.star/l

se.star <- as.numeric(sqrt(crossprod(gradient, Psi.hat.star)%*%gradient/T))
d.star  <- abs(D.hat.star - D.hat)/se.star

if(d.star >= d){p.value <- p.value + 1}
}
cat("\n ========================== \n ")

p.value = p.value/(M + 1)

out <- list("Diff"=D.hat, "p.value"=p.value)
return(out)
}

##########################
se.Parzen <- function(ret, est) 
{
if(!(est %in% c("var", "sr"))){stop("\n est must be either 'var' or 'sr'\n\n")} else
if(est=="sr") {T <- NROW(ret); obj <- sr.util(ret)} else
if(est=="var"){T <- NROW(ret); obj <- var.util(ret)}

gradient <- obj$grad; V.hat <- obj$V.hat
Psi.hat  <- Psi.hat.fn(V.hat)

return(as.numeric(sqrt(crossprod(gradient, Psi.hat)%*%gradient/T) ) )
}

##########################
se.Parzen.pw <- function(ret, est) 
{
if(!(est %in% c("var", "sr"))){stop("\n est must be either 'var' or 'sr'\n\n")} else
if(est=="sr") {T <- NROW(ret); obj <- sr.util(ret)} else
if(est=="var"){T <- NROW(ret); obj <- var.util(ret)}

gradient <- obj$grad; V.hat <- obj$V.hat;
tmp <- prewhite.fn(V.hat); V.star <- tmp$V.star; D <- tmp$D

Psi.hat <- Psi.hat.fn(V.star); Psi.hat <- D%*%tcrossprod(Psi.hat, D)

return(as.numeric(sqrt(crossprod(gradient, Psi.hat)%*%gradient/T)))
}

##########################
prewhite.fn <- function(V.hat)
{

T <- NROW(V.hat); seq1 <- seq(1, T-1); seq2 <- seq(2, T)

A.ls <- matrix(0, 4, 4); V.star <- matrix(0, T-1, 4)
reg1 <- V.hat[seq1, 1]; reg2 <- V.hat[seq1, 2]
reg3 <- V.hat[seq1, 3]; reg4 <- V.hat[seq1, 4]

# VAR FIT
for(j in (1:4))
{
fit      <- lm(V.hat[seq2,j] ~ -1 + reg1 + reg2 + reg3 + reg4)
A.ls[j,] <- as.numeric(fit$coef); V.star[, j] <- as.numeric(fit$resid)
}

svd.A <- svd(A.ls); d <- svd.A$d; d.adj <- d

# VAR COEFS ADJ
for(i in (1:4))
{
if(d[i] >  0.97){d.adj[i] <- 0.97} else
if(d[i] < -0.97){d.adj[i] <- -0.97}
}

A.hat   <- svd.A$u%*%tcrossprod(diag(d.adj), svd.A$v)
D       <- solve(diag(4) - A.hat)
reg.mat <- rbind(reg1, reg2, reg3, reg4)

for(j in (1:4)){V.star[,j] <- V.hat[seq2,j] - A.hat[j,]%*%reg.mat}

out <- list("V.star"=V.star, "D"=D)
return(out)
}

##########################
Psi.hat.fn <- function(V.hat) 
{
T <- NROW(V.hat)

alpha.hat <- alpha.hat.fn(V.hat)
S.star    <- 2.6614*(alpha.hat*T)^0.2
S.star    <- min(S.star, T-1)
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

#####################################
sr.util <- function(ret)
{
# utility function to calculate the V.hat AND the gradient
dat <- cbind(ret, ret^2); mu <- apply(dat, 2, mean)

gradient <- rep(NA, 4)
gradient[1] <-  mu[3]/(mu[3]-mu[1]^2)^1.5
gradient[2] <- -mu[4]/(mu[4]-mu[2]^2)^1.5
gradient[3] <- -(0.5)*mu[1]/(mu[3]-mu[1]^2)^1.5
gradient[4] <-  (0.5)*mu[2]/(mu[4]-mu[2]^2)^1.5

V.hat <- cbind(dat[,1]-mu[1], dat[,2]-mu[2], dat[,3]-mu[3], dat[,4]-mu[4])

out <- list("grad"=gradient, "V.hat"=V.hat)
return(out)
}

#####################################
var.util <- function(ret)
{
# utility function to calculate the V.hat AND the gradient
dat <- cbind(ret, ret^2); mu <- apply(dat, 2, mean)

gradient <- rep(NA, 4)
gradient[1] <- -2*mu[1]/(mu[3] - mu[1]^2)
gradient[2] <-  2*mu[2]/(mu[4] - mu[2]^2)
gradient[3] <-  1/(mu[3] - mu[1]^2)
gradient[4] <- -1/(mu[4] - mu[2]^2)

V.hat <- cbind(dat[,1]-mu[1], dat[,2]-mu[2], dat[,3]-mu[3], dat[,4]-mu[4])

out <- list("grad"=gradient, "V.hat"=V.hat)
return(out)
}

#####################################
cbb.seq <- function(T, b)
{
# Circluar bootstrap # Boot data is made of l blocks of size b
# circular bootstrap has UNIT MASS (only b=5 allowed)
T1 <- seq(1, T); l <- floor(T/b);

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
# sb.seq (MEU)
sb.seq <- function(T, b.av) 
{
# Stationary bootstrap # Boot data is made of l block with average size b

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

#####################################
