#######################################
# Mean Testing

#######################################
# Data and Functions
ret <- readRDS("./DATA/rets.RDS")
source("./R/LW-SR-VAR.R")

#######################################
# HEDGE FUND DATASET

# Built in Function
t.test(ret$hedge[,1], ret$hedge[,2], paired=T, alternative = "two.sided")

# My t-test
unlist(mu.t.test(ret$hedge))

# HAC test
unlist(mu.hac.test(ret$hedge))

# Boot test
mu.boot.test(ret$hedge)

#######################################
# MUTAL FUND DATASET

# Built in Function
t.test(ret$mutual[,1], ret$mutual[,2], paired=T, alternative = "two.sided")

# My t-test
unlist(mu.t.test(ret$mutual))

# HAC test
unlist(mu.hac.test(ret$mutual))

# boot test
mu.boot.test(ret$mutual)

#######################################
cat(" \n ***** END OF SCRIPT ****** \n")
