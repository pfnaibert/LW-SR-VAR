Replicating Ledoit and Wolf (2008, 2011) Functions
================
Paulo Ferreira Naibert
14, outubro, 2019

Introduction
============

On Ledoit and Wolf (2008) and Ledoit and Wolf (2011), the authors provide a method to test the hypothesis of difference in the Sharpe Ratio and the Variance of two investment strategies. The papers can be found in the following links:

[SR Paper](https://www.econ.uzh.ch/dam/jcr:ffffffff-935a-b0d6-0000-00007214c2bc/jef_2008pdf.pdf), [SR WP](http://www.econ.uzh.ch/static/wp_iew/iewwp320.pdf), [Variance Paper](https://www.econ.uzh.ch/dam/jcr:520edf26-2322-4708-8dde-51d61141914a/Ledoit_et_al-2011-Wilmott_Robust_Performance.pdf), [Variance WP](http://www.econ.uzh.ch/static/wp_iew/iewwp516.pdf).

In this post we are going to replicate Ledoit and Wolf (2008) and Ledoit and Wolf (2011). The codes can be found in [Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html), in the [publications section](https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html). Our goal is to verify if I replicated the functions on Michael Wolf's website. On what follows is a test of the functions' outputs.

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
