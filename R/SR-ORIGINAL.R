##########################
# ORIGINAL LW CODE for SHARPE RATIO
# Minor Modifications
##########################

##########################
boot.time.inference <- function (ret, b=5, M=999, Delta.null = 0, digits = 4) 
{
T <- NROW(ret); l <- floor(T/b)
Delta.hat <- sharpe.ratio.diff(ret)
d <- abs(Delta.hat - Delta.null)/compute.se.Parzen.pw(ret)
p.value <- 1

for(m in (1:M)){
ret.star <- ret[cbb.sequence(T, b),]; 
Delta.hat.star <- sharpe.ratio.diff(ret.star)

ret1.star       <- ret.star[, 1];     ret2.star       <- ret.star[, 2]
mu1.hat.star    <- mean(ret1.star);   mu2.hat.star    <- mean(ret2.star)
gamma1.hat.star <- mean(ret1.star^2); gamma2.hat.star <- mean(ret2.star^2)

gradient <- rep(0, 4)
gradient[1] <-  gamma1.hat.star/(gamma1.hat.star - mu1.hat.star^2)^1.5
gradient[2] <- -gamma2.hat.star/(gamma2.hat.star - mu2.hat.star^2)^1.5
gradient[3] <- -0.5 * mu1.hat.star/(gamma1.hat.star - mu1.hat.star^2)^1.5
gradient[4] <-  0.5 * mu2.hat.star/(gamma2.hat.star - mu2.hat.star^2)^1.5

y.star <- data.frame(ret1.star - mu1.hat.star, ret2.star - mu2.hat.star, ret1.star^2 - gamma1.hat.star, ret2.star^2 - gamma2.hat.star)

Psi.hat.star <- matrix(0, 4, 4)
for(j in (1:l))
{
zeta.star    <- b^0.5 * colMeans(y.star[(1 + (j - 1)*b):(j*b), ])
Psi.hat.star <- Psi.hat.star + zeta.star%*%t(zeta.star)
}
Psi.hat.star <- Psi.hat.star/l

se.star <- as.numeric(sqrt(t(gradient)%*%Psi.hat.star%*%gradient/T))
d.star  <- abs(Delta.hat.star - Delta.hat)/se.star

if(d.star >= d){p.value <- p.value + 1}
}

p.value <- p.value/(M + 1)

out <- list(Difference = round(Delta.hat, digits), p.Value = round(p.value,digits))
return(out)
}

##########################
sharpe.ratio.diff <- function(ret) 
{
ret1     <- ret[, 1];   ret2     <- ret[, 2]
mu1.hat  <- mean(ret1); mu2.hat  <- mean(ret2)
sig1.hat <- sd(ret1);   sig2.hat <- sd(ret2)

diff = mu1.hat/sig1.hat - mu2.hat/sig2.hat

return(diff)
}

##########################
hac.inference <- function(ret, digits = 3) 
{
ret1 = ret[, 1]
ret2 = ret[, 2]

mu1.hat = mean(ret1)
mu2.hat = mean(ret2)

sig1.hat = sd(ret1)
sig2.hat = sd(ret2)

SR1.hat = mu1.hat / sig1.hat
SR2.hat = mu2.hat / sig2.hat

SRs = c(SR1.hat, SR2.hat); names(SRs) = c("SR1.hat", "SR2.hat")

diff = SR1.hat - SR2.hat

se = compute.se.Parzen(ret)
se.pw = compute.se.Parzen.pw(ret)
SEs = c(se, se.pw); names(SEs) = c("HAC", "HAC.pw")

PV = 2 * pnorm(-abs(diff)/se)
PV.pw = 2 * pnorm(-abs(diff)/se.pw)
PVs = round(c(PV, PV.pw), digits)
names(PVs) = c("HAC", "HAC.pw")

out <- list(Sharpe.Ratios = SRs, Difference = diff, Standard.Errors = SEs, p.Values = PVs)

return(out)
}

##########################
compute.se.Parzen <- function(ret) 
{
T = NROW(ret)

ret1 = ret[, 1]; ret2 = ret[, 2]
mu1.hat = mean(ret1); mu2.hat = mean(ret2)
ret1.2 = ret1^2; ret2.2 = ret2^2
gamma1.hat = mean(ret1.2); gamma2.hat = mean(ret2.2)

gradient = rep(0, 4)
gradient[1] = gamma1.hat/(gamma1.hat - mu1.hat^2)^1.5
gradient[2] = -gamma2.hat/(gamma2.hat - mu2.hat^2)^1.5
gradient[3] = -0.5 * mu1.hat/(gamma1.hat - mu1.hat^2)^1.5
gradient[4] = 0.5 * mu2.hat/(gamma2.hat - mu2.hat^2)^1.5

V.hat = cbind(ret1 - mu1.hat, ret2 - mu2.hat, ret1.2 - gamma1.hat, ret2.2 - gamma2.hat)
Psi.hat = compute.Psi.hat(V.hat)

se = as.numeric(sqrt(t(gradient) %*% Psi.hat %*% gradient/T))
        
return(se)
}

##########################
compute.se.Parzen.pw <- function (ret) 
{
T = NROW(ret)

ret1 = ret[, 1]; ret2 = ret[, 2]
mu1.hat = mean(ret1); mu2.hat = mean(ret2)
ret1.2 = ret1^2; ret2.2 = ret2^2
gamma1.hat = mean(ret1.2); gamma2.hat = mean(ret2.2)

gradient = rep(0, 4)
gradient[1] = gamma1.hat/(gamma1.hat - mu1.hat^2)^1.5
gradient[2] = -gamma2.hat/(gamma2.hat - mu2.hat^2)^1.5
gradient[3] = -0.5 * mu1.hat/(gamma1.hat - mu1.hat^2)^1.5
gradient[4] = 0.5 * mu2.hat/(gamma2.hat - mu2.hat^2)^1.5


V.hat = cbind(ret1 - mu1.hat, ret2 - mu2.hat, ret1.2 - gamma1.hat, ret2.2 - gamma2.hat)
A.ls = matrix(0, 4, 4)
V.star = matrix(0, T - 1, 4)

reg1 = V.hat[1:T - 1, 1]
reg2 = V.hat[1:T - 1, 2]
reg3 = V.hat[1:T - 1, 3]
reg4 = V.hat[1:T - 1, 4]

for (j in (1:4)) {
    fit = lm(V.hat[2:T, j] ~ -1 + reg1 + reg2 + reg3 + reg4)
    A.ls[j, ] = as.numeric(fit$coef)
    V.star[, j] = as.numeric(fit$resid)
}

svd.A = svd(A.ls)
d = svd.A$d
d.adj = d

for (i in (1:4)) {
    if (d[i] > 0.97) 
        d.adj[i] = 0.97
    else if (d[i] < -0.97) 
        d.adj[i] = -0.97
}

A.hat = svd.A$u %*% diag(d.adj) %*% t(svd.A$v)
D = solve(diag(4) - A.hat)
reg.mat = rbind(reg1, reg2, reg3, reg4)

for (j in (1:4)) {
    V.star[, j] = V.hat[2:T, j] - A.hat[j, ] %*% reg.mat
}

Psi.hat = compute.Psi.hat(V.star)
Psi.hat = D %*% Psi.hat %*% t(D)
se = as.numeric(sqrt(gradient %*% Psi.hat %*% gradient/T))

return(se)
}

##########################
compute.V.hat <- function(ret) 
{
ret1 = ret[, 1]; ret2 = ret[, 2]
V.hat = cbind(ret1 - mean(ret1), ret2 - mean(ret2), ret1^2 - mean(ret1^2), ret2^2 - mean(ret2^2))
return(V.hat)
}

##########################
compute.Psi.hat <- function(V.hat) 
{
T = length(V.hat[, 1])

alpha.hat = compute.alpha.hat(V.hat)
S.star = 2.6614 * (alpha.hat * T)^0.2
S.star = min(S.star, T-1)
Psi.hat = compute.Gamma.hat(V.hat, 0)

j = 1
while (j < S.star) {
    Gamma.hat = compute.Gamma.hat(V.hat, j)
    Psi.hat = Psi.hat + kernel.Parzen(j/S.star) * (Gamma.hat + t(Gamma.hat))
    j = j + 1
}

Psi.hat = (T/(T - 4))*Psi.hat

return(Psi.hat)
}

##########################
compute.alpha.hat <- function(V.hat) 
{

T <- NROW(V.hat); p <- NCOL(V.hat)
num <- 0; den <- 0

for (i in (1:p))
{
fit = ar(V.hat[, i], 0, 1, method = "ols")
rho.hat = as.numeric(fit[2])
sig.hat = sqrt(as.numeric(fit[3]))
num = num + 4 * rho.hat^2 * sig.hat^4/(1 - rho.hat)^8
den = den + sig.hat^4/(1 - rho.hat)^4
}

return(num/den)
}

##########################
compute.Gamma.hat <- function (V.hat, j) 
{
T <- NROW(V.hat); p <- NCOL(V.hat); Gamma.hat <- matrix(0, p, p)

if (j >= T){ stop("j must be smaller than the row dimension!") }

for (i in ((j + 1):T)){Gamma.hat = Gamma.hat + V.hat[i, ]%*%t(V.hat[i - j, ])} 

return(Gamma.hat/T)
}

##########################
kernel.Parzen <- function(x) 
{
if (abs(x) <= 0.5){result = 1 - 6*x^2 + 6*abs(x)^3}

else if (abs(x) <= 1){result = 2*(1 - abs(x))^3}

else result = 0

return(result)
}

##########################
cbb.sequence <- function (T, b)
{
l = floor(T/b)
index.sequence = c(1:T, 1:b)
sequence = rep(0, T)
start.points = sample(1:T, l, replace = T)

for (j in (1:l)) {
    start = start.points[j]
    sequence[((j - 1) * b + 1):(j * b)] = index.sequence[start:(start + b - 1)]
}

return(sequence)
}

##########################
sb.sequence <- function(T, b.av, length = T) 
{
index.sequence = c(1:T, 1:T); sequence = rep(0, length + T)

current = 0
while (current < length) {
    start = sample(1:T, 1)
    b = rgeom(1, 1/b.av) + 1
    sequence[(current + 1):(current + b)] = index.sequence[start:(start + b - 1)]
    current = current + b
}

return(sequence[1:length])
}

##########################
block.size.calibrate <- function (ret, b.vec = c(1, 3, 6, 10), alpha = 0.05, M = 199, K = 1000, b.av = 5, T.start = 50) 
{

T = NROW(ret1); b.len = length(b.vec); emp.reject.probs = rep(0, b.len)

Delta.hat = sharpe.ratio.diff(ret)
ret1 = ret[, 1]
ret2 = ret[, 2]

Var.data = matrix(0, T.start + T, 2)
Var.data[1, 1] = ret[1, 1]
Var.data[1, 2] = ret[1, 2]

Delta.hat = sharpe.ratio.diff(ret)

fit1 = lm(ret1[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])
fit2 = lm(ret2[2:T] ~ ret1[1:(T - 1)] + ret2[1:(T - 1)])

coef1 = as.numeric(fit1$coef)
coef2 = as.numeric(fit2$coef)

resid.mat = cbind(as.numeric(fit1$resid), as.numeric(fit2$resid))

for (k in (1:K))
{
resid.mat.star = rbind(c(0, 0), resid.mat[sb.sequence(T - 1, b.av, T.start + T - 1), ])
# print(resid.mat.star)

for (t in (2:(T.start + T)))
{
Var.data[t, 1] = coef1[1] + coef1[2] * Var.data[t - 1, 1] + coef1[3] * Var.data[t - 1, 2] + resid.mat.star[t,1]

Var.data[t, 2] = coef2[1] + coef2[2] * Var.data[t - 1, 1] + coef2[3] * Var.data[t - 1, 2] + resid.mat.star[t,2]
}
# print(Var.data)

Var.data.trunc = Var.data[(T.start + 1):(T.start + T),]

for (j in (1:b.len))
{
p.Value = boot.time.inference(Var.data.trunc, b.vec[j], M, Delta.hat)$p.Value

if (p.Value <= alpha){emp.reject.probs[j] = emp.reject.probs[j] + 1 }
}

}

emp.reject.probs = emp.reject.probs/K
b.order = order(abs(emp.reject.probs - alpha))
b.opt = b.vec[b.order[1]]
b.vec.with.probs = rbind(b.vec, emp.reject.probs)
colnames(b.vec.with.probs) = rep("", length(b.vec))

out <- list(Empirical.Rejection.Probs = b.vec.with.probs, b.optimal = b.opt)

return(out)
}

#####################################
