Replicating Ledoit and Wolf (2008, 2011) Functions
================
Paulo Ferreira Naibert

Original: 2019 October 13. Current: 14, outubro, 2019

Introduction
============

On Ledoit and Wolf (2008) and Ledoit and Wolf (2011), the authors provide a method to test the hypothesis of difference in the Sharpe Ratio and the Variance of two investment strategies. The papers can be found in the following links:

[SR Paper](https://www.econ.uzh.ch/dam/jcr:ffffffff-935a-b0d6-0000-00007214c2bc/jef_2008pdf.pdf), [SR WP](http://www.econ.uzh.ch/static/wp_iew/iewwp320.pdf), [Variance Paper](https://www.econ.uzh.ch/dam/jcr:520edf26-2322-4708-8dde-51d61141914a/Ledoit_et_al-2011-Wilmott_Robust_Performance.pdf), [Variance WP](http://www.econ.uzh.ch/static/wp_iew/iewwp516.pdf).

In this post we are going to replicate Ledoit and Wolf (2008) and Ledoit and Wolf (2011). The codes can be found in [Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html), in the [publications section](https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html). Our goal is to verify if I replicated the functions on Michael Wolf's website. On what follows is a test of the functions' outputs.

Functions
---------

I copied the original functions on `Sharpe.RData` and `Var.RData` into the files: `"./SR-ORIGINAL.R"` and `"./Var-ORIGINAL.R"` with some minor modifications. Then I copied the functions to the file `"./MYfuns.R"` and modified them more.

Data
----

For the data, I simpy saved the returns data of `Sharpe.RData` as a `.RDS` file. The following code should work:

``` r
load("Sharpe.RData")
lwret <- list("mutual"=as.matrix(ret.agg), "hedge"=ret.hedge)
saveRDS(lwret, "./rets.RDS")
```

Loading `rets.RDS` should result in a list object on the `R` Software (R Core Team 2019) .

Functions on Each File
======================

First let's take a look on which funs the files share.

``` r
# Differences in the funs

# Funs in "./SR-ORIGINAL.R"
source("./SR-ORIGINAL.R");
sr <- ls();
rm(list=setdiff(ls(), "sr"))
sr
```

    ##  [1] "block.size.calibrate" "boot.time.inference"  "cbb.sequence"        
    ##  [4] "compute.alpha.hat"    "compute.Gamma.hat"    "compute.Psi.hat"     
    ##  [7] "compute.se.Parzen"    "compute.se.Parzen.pw" "compute.V.hat"       
    ## [10] "hac.inference"        "kernel.Parzen"        "sb.sequence"         
    ## [13] "sharpe.ratio.diff"

``` r
# Funs in "./Var-ORIGINAL.R"
source("./Var-ORIGINAL.R");
lv <- setdiff(ls(), "sr");
rm(list=setdiff(ls(), c("sr", "lv"))); 
lv
```

    ##  [1] "block.size.calibrate.log.var" "boot.time.inference.log.var" 
    ##  [3] "cbb.sequence"                 "compute.alpha.hat"           
    ##  [5] "compute.Gamma.hat"            "compute.Psi.hat"             
    ##  [7] "compute.se.Parzen.log.var"    "compute.se.Parzen.pw.log.var"
    ##  [9] "compute.V.hat"                "hac.inference.log.var"       
    ## [11] "kernel.Parzen"                "log.var.diff"                
    ## [13] "sb.sequence"

``` r
# Shared functions on files
intersect(sr, lv)
```

    ## [1] "cbb.sequence"      "compute.alpha.hat" "compute.Gamma.hat"
    ## [4] "compute.Psi.hat"   "compute.V.hat"     "kernel.Parzen"    
    ## [7] "sb.sequence"

``` r
# Function on SR but not in Var
setdiff(sr, lv)
```

    ## [1] "block.size.calibrate" "boot.time.inference"  "compute.se.Parzen"   
    ## [4] "compute.se.Parzen.pw" "hac.inference"        "sharpe.ratio.diff"

``` r
# Function on Var but not in SR
setdiff(lv, sr)
```

    ## [1] "block.size.calibrate.log.var" "boot.time.inference.log.var" 
    ## [3] "compute.se.Parzen.log.var"    "compute.se.Parzen.pw.log.var"
    ## [5] "hac.inference.log.var"        "log.var.diff"

``` r
# Funs in "./MYfuns.R"
source("./MYfuns.R"); my <- setdiff(ls(), c("sr","lv") )
my
```

    ##  [1] "alpha.hat.fn"  "boot.inf.sr"   "boot.inf.var"  "cbb.seq"      
    ##  [5] "Gamma.hat.fn"  "hac.inf.sr"    "hac.inf.var"   "kernel.Parzen"
    ##  [9] "lv.diff"       "prewhite.fn"   "Psi.hat.fn"    "sb.seq"       
    ## [13] "se.Parzen"     "se.Parzen.pw"  "sr.diff"       "sr.util"      
    ## [17] "var.util"

``` r
# Shared functions on files
intersect(my, sr)
```

    ## [1] "kernel.Parzen"

``` r
intersect(my, lv)
```

    ## [1] "kernel.Parzen"

``` r
# remove all funs
rm(list=ls())
```

So I just didn't change the name of `kernel.Parzen()`

Data and Functions
==================

First, let's check if there are not any loaded objects

``` r
ls()
```

    ## character(0)

Now we load the data and check it.

``` r
# Load Data
ret <- readRDS("./rets.RDS")
str(ret)
```

    ## List of 2
    ##  $ mutual: num [1:120, 1:2] 3.93 -2.33 -4.86 1.96 -0.44 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:120] "1" "2" "3" "4" ...
    ##   .. ..$ : chr [1:2] "V1" "V2"
    ##  $ hedge : num [1:120, 1:2] 0.681 0.659 1.113 0.525 0.825 ...

Only now do we load the functions.

``` r
# Load my functions
source("./MYfuns.R")

# Load LW functions
source("./SR-ORIGINAL.R"); source("./Var-ORIGINAL.R");
```

HAC INFERENCE
=============

Here, and in the next sections, we use the `all.equal()` function with argument `check.names` set to `FALSE`. Now let's look if the output of the functions are different.

Sharpe Ratio
------------

Starting with the Sharpe Ratio:

``` r
tmp1 <- hac.inf.sr(ret$hedge, est="sr")
tmp2 <- hac.inference(ret$hedge)

tmp1; tmp2; all.equal(tmp1, tmp2, check.names=FALSE)
```

    ## $SRs
    ##      SR1      SR2 
    ## 1.014228 1.460547 
    ## 
    ## $Diff
    ## [1] -0.4463183
    ## 
    ## $SEs
    ##       HAC    HAC.pw 
    ## 0.3176271 0.3872694 
    ## 
    ## $pVals
    ##       HAC    HAC.pw 
    ## 0.1599724 0.2491259

    ## $Sharpe.Ratios
    ##  SR1.hat  SR2.hat 
    ## 1.014228 1.460547 
    ## 
    ## $Difference
    ## [1] -0.4463183
    ## 
    ## $Standard.Errors
    ##       HAC    HAC.pw 
    ## 0.3176271 0.3872694 
    ## 
    ## $p.Values
    ##    HAC HAC.pw 
    ##  0.160  0.249

    ## [1] "Component 4: Mean relative difference: 0.0003750761"

The only difference is on component 4. This is the \(p\)-Values object. The difference is on the 4th decimal place. The difference might be due to the `round()` function. On `SR-ORIGINAL.R` this value is set with `round(,3)`

Variances
---------

Now, the variances:

``` r
tmp1 <- hac.inf.var(ret$hedge, est="var")
tmp2 <- hac.inference.log.var(ret$hedge)

tmp1; tmp2; all.equal(tmp1, tmp2, check.names=FALSE)
```

    ## $Vars
    ##       Var1       Var2 
    ## 1.46589237 0.02809697 
    ## 
    ## $Log.Vars
    ##       Var1       Var2 
    ##  0.3824642 -3.5720937 
    ## 
    ## $Diff
    ## [1] 3.954558
    ## 
    ## $SEs
    ##      HAC   HAC.pw 
    ## 0.514087 0.699845 
    ## 
    ## $pVals
    ##          HAC       HAC.pw 
    ## 1.444110e-14 1.598705e-08

    ## $Variances
    ## Var1.hat Var2.hat 
    ##    1.466    0.028 
    ## 
    ## $Log.Variances
    ## LogVar1.hat LogVar2.hat 
    ##       0.382      -3.572 
    ## 
    ## $Difference
    ## [1] 3.955
    ## 
    ## $Standard.Errors
    ##    HAC HAC.pw 
    ##  0.514  0.700 
    ## 
    ## $p.Values
    ##    HAC HAC.pw 
    ##      0      0

    ## [1] "Component 1: Mean relative difference: 0.0001369429"
    ## [2] "Component 2: Mean relative difference: 0.0001410697"
    ## [3] "Component 3: Mean relative difference: 0.0001118031"
    ## [4] "Component 4: Mean relative difference: 0.0001993381"

All outputs are different: `Variances` `Log.Variances` `Difference` `Standard.Errors` and `p.Values` However, we can note the the differences are only on the 4th decimal place. This might be due to the `round()` function. On `Var-ORIGINAL.R` this value is set with `round(,3)`

BOOT TIME SERIES INFERENCE
==========================

Now we turn to the Inference using Bootstrap. We first set a seed. The reasons for doing so is twofold: (a) so we can repeat our experiment, and (b) so the output of the functions are comparable.

``` r
# Set seed so results are comparable
set.seed(1)
```

Then, we evaluate the difference in the Sharpe Ratio, again using the `all.equal()` function.

Sharpe Ratio
------------

``` r
tmp1 <- boot.inf.sr(ret$hedge, est="sr", b=5, M=499)
tmp2 <- boot.time.inference(ret$hedge, b=5, M=499)

tmp1; tmp2; all.equal(tmp1, tmp2, check.names=FALSE)
```

    ## $Diff
    ## [1] -0.4463
    ## 
    ## $pVal
    ## [1] 0.298

    ## $Difference
    ## [1] -0.4463
    ## 
    ## $p.Value
    ## [1] 0.25

    ## [1] "Component 2: Mean relative difference: 0.1610738"

We can note that we have a difference in the p-values of the Sharpe Ratios. This might be due to the small sample correction \(T/(T-4)\) on the `Psi.Star` quantity that I introduced in my code, but not on LW's.

Variance
--------

Finally, we check the output of the function that performs boostrap inference in the variance.

``` r
tmp1 <- boot.inf.var(ret$hedge, est="var", b=5, M=499)
tmp2 <- boot.time.inference.log.var(ret$hedge, b=5, M=499)

tmp1; tmp2; all.equal(tmp1, tmp2, check.names=FALSE)
```

    ## $Difference
    ## [1] 3.955
    ## 
    ## $p.Value
    ## [1] 0.002

    ## $Difference
    ## [1] 3.955
    ## 
    ## $p.Value
    ## [1] 0.002

    ## [1] TRUE

The output of the `all.equal()` function is `TRUE` which indicates that all quantities are equal.

System Information
==================

Software Information
--------------------

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=pt_BR.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=pt_BR.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=pt_BR.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=pt_BR.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] tufterhandout_1.2.1 tufte_0.5           rmarkdown_1.16     
    ## [4] knitr_1.25         
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.6.1  magrittr_1.5    htmltools_0.4.0 tools_3.6.1    
    ##  [5] yaml_2.2.0      Rcpp_1.0.2      stringi_1.4.3   stringr_1.2.0  
    ##  [9] digest_0.6.21   xfun_0.10       rlang_0.4.0     evaluate_0.14

DISCLAIMER
==========

I own NONE of the rights to the codes. I simply edited the already existing codes in [Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html). I also offer NO support for the functions. They are presented AS IS.

REFERENCES
==========

Ledoit, Oliver, and Michael Wolf. 2008. “Robust performance hypothesis testing with the Sharpe ratio.” *Journal of Empirical Finance* 15 (5): 850–59.

———. 2011. “Robust performance hypothesis testing with the Variance.” *Wilmott Magazine*, no. 55 (September): 86–89.

R Core Team. 2019. *R: A Language and Environment for Statistical Computing*. Vienna, Austria: R Foundation for Statistical Computing. <https://www.R-project.org/>.
