#######################################
library(knitr)
ret <- readRDS("./DATA/rets.RDS")

#######################################
# TABS
d.h <- ret$hedge[,1]-ret$hedge[,2]
d.m <- ret$mutual[,1]-ret$mutual[,2]

mus <- c(mean(ret$hedge[,1]), mean(ret$hedge[,2]), mean(d.h),
        mean(ret$mutual[,1]), mean(ret$mutual[,2]), mean(d.m) )
sds <- c(sd(ret$hedge[,1]), sd(ret$hedge[,2]), sd(d.h),
        sd(ret$mutual[,1]), sd(ret$mutual[,2]), sd(d.m) )
srs <- mus/sds

tab <- cbind(mus, sds, srs); rownames(tab) <- c("hedge1", "hedge2", "diff.hedge", "mutual1", "mutual2", "diff.mutual")
kable(tab, digits=3)

#######################################
# PLOTS

par(mfrow=c(3,1))
plot(density(ret$hedge[,1]), main="Density Plot of hedge1")
plot(density(ret$hedge[,2]), main="Density Plot of hedge2")
plot(density(d.h), main="Density Plot of Ret Diff of Hedges")

par(mfrow=c(3,1))
plot(density(ret$mutual[,1]), main="Density Plot of mutual1")
plot(density(ret$mutual[,2]), main="Density Plot of mutual2")
plot(density(d.m), main="Density Plot of Ret Diff of Mutuals")

#######################################
# TODO
# rets against time
# rets^2 against time
# autocorrelation

#######################################
cat(" \n ***** END OF SCRIPT ****** \n")
