#' --- 
#' title: DEMO of Replicating Ledoit and Wolf (2008, 2011) Functions
#' author: "Paulo Ferreira Naibert" 
#' output: 
#'   html_document:
#'     toc: false
#'     number_sections: true
#' keep_tex: true
#' bibliography: refs.bib
#' --- 

# /* ################################################# */
#+ setup, include=FALSE
knitr::opts_chunk$set(error = TRUE)
# remove all objects
rm(list=ls()) 

# /* ################################################# */
#'Original: 2019 October 13.
#'
#'Current: `r format(Sys.time(), "%d, %B, %Y")`
#'
#' # Overview
#' The goal of this DEMO to verify if the repository replicates the functions on Michael Wolf's website.
#' 
#' I didn't have time to replicate the `block.size.calibrate()` functions.
#' 
#' The ORIGINAL codes can be found in 
#' [Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html), 
#' in the 
#' [publications section](https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html).
#'

# /* ################################################# */
#' # Differences in the Functions on Each File
#'
#' First let's take a look on which funs the files share.
#' 

# /* ################################################# */
#' ## Funs in `Shape.RData`
ls()
load("./R/Sharpe.RData");
sr <- ls();
rm(list=setdiff(ls(), "sr"))
sr

# /* ################################################# */
#' ## Funs in `Var.RData`
load("./R/Var.RData");
lv <- setdiff(ls(), "sr");
rm(list=setdiff(ls(), c("sr", "lv"))); 
lv

#' ### Shared functions on files
intersect(sr, lv)

#' ### Function on SR but not in Var
setdiff(sr, lv)

#' ### Function on Var but not in SR
setdiff(lv, sr)

# /* ################################################# */
#' ## Funs in `./R/LW-SR.Var.R`
source("./R/LW-SR-VAR.R"); my <- setdiff(ls(), c("sr","lv") )
my

#' ### Shared functions on files
intersect(my, sr)
intersect(my, lv)

#' So I just didn't change the name of `kernel.Parzen()`.
#' Let's check the usage of the functions.
#'

# /* ################################################# */
#' # Data and Functions
#'
#' First, we remove all objects and then we  check if there are not any loaded objects
rm(list=ls())
ls()

#' Now we load the data and check it.
ret <- readRDS("./DATA/rets.RDS")
str(ret)

#' Only now do we load the functions.
#'
source("./R/LW-SR-VAR.R");

#' Load LW functions
#'
load("./R/Var.RData"); load("./R/Sharpe.RData");

# /* ################################################# */
#' est.diff
est.diff(ret$hedge, est="mu")
est.diff(ret$hedge, est="var")
est.diff(ret$hedge, est="sr")
est.diff(ret$hedge, est="ceq")
est.diff(ret$hedge, est="ceq", rav=1)
est.diff(ret$hedge, est="ceq", rav=2)
est.diff(ret$hedge, est="foo")

# /* ################################################# */
#' util.fn
str( util.fn(ret$hedge, est="mu")  )
str( util.fn(ret$hedge, est="var") )
str( util.fn(ret$hedge, est="sr")  )
str( util.fn(ret$hedge, est="ceq") )
str( util.fn(ret$hedge, est="foo") )

# /* ################################################# */
#' # Classic Testing
#' `diff.test(ret, est, rav=1, kernel="parzen", pw=0)` function.
#' 
#' In classic testing, we use `kernel="sample":
#' 
#' Classic Tests:
diff.test(ret$hedge, est = "mu", kernel="sample")
#' Let's compare with the `t.test` built-in function:
t.test(ret$hedge[,1], ret$hedge[,2], paired=T, alternative = "two.sided")
#' Comparing the test value of both tests, we can see that they are pretty close.
#' 
#' Other Statistics:
diff.test(ret$hedge, est = "var", kernel="sample")
diff.test(ret$hedge, est = "sr", kernel="sample")
diff.test(ret$hedge, est = "ceq", kernel="sample")
diff.test(ret$hedge, est = "ceq", rav=3, kernel="sample")

#' Let's see the error messages the functions return:
#' ERROR MESSAGES
diff.test(ret$hedge, est = "foo", rav=3, kernel="sample")
diff.test(ret$hedge, est = "ceq", rav=3, kernel="foo")
diff.test(ret$hedge, est = "foo", rav=3, kernel="foo")

# /* ################################################# */
#' #  HAC INFERENCE
#' Here, and in the next sections, we use the `all.equal()` function with argument `check.names` set to `FALSE` to attest if our outputs diverge from the LW functions.
#'
#' ## Sharpe Ratio
#' `diff.test(ret, est, rav=1, kernel="parzen", pw=0)` function.
tmp1 <- unlist(hac.inference(ret$hedge)); tmp1
tmp2 <- unlist(diff.test(ret$hedge, est="sr", kernel="parzen", pw=0)); tmp2
tmp3 <- unlist(diff.test(ret$hedge, est="sr", kernel="parzen", pw=1)); tmp3

#' ### Test SR HAC
#' Difference in `Difference`
all.equal(tmp1["Difference"], tmp2["Diff"], check.names=FALSE)
all.equal(tmp1["Difference"], tmp3["Diff"], check.names=FALSE)
#' The differences are on the 4th decimal place.
#' The difference might be due to the `round()` function.
#' On the original functions this value is set with `round(,3)`.
#'

#' Difference in `Standard.Errors`
all.equal(tmp1["Standard.Errors.HAC"], tmp2["se"], check.names=FALSE)
all.equal(tmp1["Standard.Errors.HAC.pw"], tmp3["se"], check.names=FALSE)
#' The differences are on the 3th and 4th decimal place.
#' 

#' Difference in `p.Values`
all.equal(tmp1["p.Values.HAC"], tmp2["p.value"], check.names=FALSE)
all.equal(tmp1["p.Values.HAC.pw"], tmp3["p.value"], check.names=FALSE)
#' The differences are on the 4th decimal place.


# /* ################################################# */
#' ## Variances
#' `diff.test(ret, est, rav=1, kernel="parzen", pw=0)` function.
#' Now, the variances:
tmp1 <- unlist(hac.inference.log.var(ret$hedge)); tmp1
tmp2 <- unlist(diff.test(ret$hedge, est="var", kernel="parzen", pw=0)); tmp2
tmp3 <- unlist(diff.test(ret$hedge, est="var", kernel="parzen", pw=1)); tmp3

#' ### Test Var HAC
#' Difference in `Difference`
all.equal(tmp1["Difference"], tmp2["Diff"], check.names=FALSE)
all.equal(tmp1["Difference"], tmp3["Diff"], check.names=FALSE)
#' The differences are on the 4th decimal place.
#' The difference might be due to the `round()` function.
#' On the original functions this value is set with `round(,3)`.
#'

#' Difference in `Standard.Errors`
all.equal(tmp1["Standard.Errors.HAC"], tmp2["se"], check.names=FALSE)
all.equal(tmp1["Standard.Errors.HAC.pw"], tmp3["se"], check.names=FALSE)
#' The differences are on the 4th decimal place.
#' 

#' Difference in `p.Values`
all.equal(tmp1["p.Values.HAC"], tmp2["p.value"], check.names=FALSE)
all.equal(tmp1["p.Values.HAC.pw"], tmp3["p.value"], check.names=FALSE)
#' The differences are smaller than the 4th decimal place.
#' 

# /* ################################################# */
#' ## Other Statistics

#' ### not PW
diff.test(ret$hedge, est = "mu")
diff.test(ret$hedge, est = "ceq")
diff.test(ret$hedge, est = "ceq", rav=3)

#' ### PW
diff.test(ret$hedge, est = "mu", pw=1)
diff.test(ret$hedge, est = "ceq", pw=1)
diff.test(ret$hedge, est = "ceq", rav=3, pw=1)

#' ### ERROR MESSAGES
diff.test(ret$hedge, est = "foo", rav=3)
diff.test(ret$hedge, est = "foo", rav=3, pw=1)


# /* ################################################# */
#' # BOOTSTRAP TIME SERIES INFERENCE

#' Now we turn to the Inference using Bootstrap.
#' We first use the `set.seed()` function.
#' The reasons for doing so is twofold:
#' (a) so we can repeat our experiment, and (b) so the output of the functions are comparable.
#'
#' Then, we evaluate the difference in the Sharpe Ratio, again using the `all.equal()` function.
#' 
#' ## Sharpe Ratio Boot
#' boot.test(ret, est, rav=1, b=5, M=499, D.null=0) 
#' Original:
set.seed(1)
t0 <- Sys.time()
tmp1 <- boot.time.inference(ret$hedge, b=5, M=499)
t1 <- Sys.time(); t1-t0
tmp1

#' Modified:
set.seed(1)
t0 <- Sys.time()
tmp2 <- boot.test(ret$hedge, est="sr",  b=5, M=499)
t1 <- Sys.time(); t1-t0
tmp2

#'
#' ### Test SR boot
all.equal(tmp1["Difference"], tmp2["Diff"], check.names=FALSE)
all.equal(tmp1["p.Value"], tmp2["p.value"], check.names=FALSE)

#' The output of the test indicates difference in the `Difference` argument only on the 5th decimal place.
#' The `p.value` is equal.

# /* ################################################# */
#'
#' ## Variance Boot
#' boot.test(ret, est, rav=1, b=5, M=499, D.null=0) 
#'Finally, we check the output of the function that performs boostrap inference in the variance.
#' Original:
set.seed(1)
t0 <- Sys.time()
tmp1 <- unlist(boot.time.inference.log.var(ret$hedge, b=5, M=499))
t1 <- Sys.time(); t1-t0
tmp1

#' Modified:
set.seed(1)
t0 <- Sys.time()
tmp2 <- unlist(boot.test(ret$hedge, est="var", b=5, M=499)); tmp2
t1 <- Sys.time(); t1-t0

#' ### Test Var boot
#     
all.equal(tmp1["Difference"], tmp2["Diff"], check.names=FALSE)
all.equal(tmp1["p.Value"], tmp2["p.value"], check.names=FALSE)
#' The output of the test indicates difference in the `Difference` argument only on the 4th decimal place.
#' The `p.value` is equal.

# /* ################################################# */
#'
#'## Other Statistics 
#' boot.test(ret, est, rav=1, b=5, M=499, D.null=0) 

set.seed(1)
boot.test(ret$hedge, est = "mu")

set.seed(1)
boot.test(ret$hedge, est = "ceq")

set.seed(1)
boot.test(ret$hedge, est = "ceq", rav=3)

#' ### ERROR MESSAGES:
set.seed(1)
boot.test(ret$hedge, est = "foo")

# /* ################################################# */
#' # Software Information
sessionInfo()

#######################################
#+ end, include=FALSE
cat(" \n ***** END OF SCRIPT ****** \n")
