Replication of Ledoit and Wolf (2008, 2011) Functions
================
Paulo Ferreira Naibert

Original: 2019 October 13.

Current: 18, outubro, 2019

Overview
========

This Repository contains functions to replicate Oliver Ledoit and Wolf (2008) and Oliver Ledoit and Wolf (2011) Olivier Ledoit and Wolf (2018).

Michael Wolf provides `R` codes for those two papers on [his webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html), in the [publications section](https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html).

What is found here is merely an edit of those codes.

If you find any bug, or technical error, please contact me.

PAPERS and GOALS of the code
============================

On Oliver Ledoit and Wolf (2008), Oliver Ledoit and Wolf (2011), and Olivier Ledoit and Wolf (2018) the authors provide methods to test the hypothesis of difference in the Sharpe Ratio and the Variance of two investment strategies. The papers can be found in the following links:

[SR Paper](https://www.econ.uzh.ch/dam/jcr:ffffffff-935a-b0d6-0000-00007214c2bc/jef_2008pdf.pdf), [SR WP](http://www.econ.uzh.ch/static/wp_iew/iewwp320.pdf), [Variance Paper](https://www.econ.uzh.ch/dam/jcr:520edf26-2322-4708-8dde-51d61141914a/Ledoit_et_al-2011-Wilmott_Robust_Performance.pdf), [Variance WP](http://www.econ.uzh.ch/static/wp_iew/iewwp516.pdf), [Generalized WP](http://www.econ.uzh.ch/static/wp/econwp305.pdf).

Functions
=========

I copied the original functions on `Sharpe.RData` and `Var.RData` into the files: `"SR-ORIGINAL.R"` and `"Var-ORIGINAL.R"` with some minor modifications. Then I copied the functions to the file `"LW-SR-Var.R"` and modified them more.

I didn't have enought time to modify the `block.size.calibrate()` functions.

Data
====

For the data, I simpy saved the returns data of `Sharpe.RData` as a `.RDS` file. Loading `rets.RDS` should result in a list object on the `R` Software (R Core Team 2019) .

USAGE
=====

See DEMO files.

TODO
====

-   Document Functions (R and html nb)

DISCLAIMER
==========

I own NONE of the rights to the codes. I simply edited the already existing codes in [Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html). I also offer NO support for the functions. They are presented "AS IS".

ACKNOWLEDGEMENTS
================

I would like to thank Michael Wolf for providing the original functions on his webstite.

REFERENCES
==========

Ledoit, Oliver, and Michael Wolf. 2008. “Robust performance hypothesis testing with the Sharpe ratio.” *Journal of Empirical Finance* 15 (5): 850–59.

———. 2011. “Robust performance hypothesis testing with the Variance.” *Wilmott Magazine*, no. 55 (September): 86–89.

Ledoit, Olivier, and Michael Wolf. 2018. “Robust Performance Hypothesis Testing with Smooth Functions of Population Moments.” Working Paper Series 305. Department of Economics, University of Zurich.

R Core Team. 2019. *R: A Language and Environment for Statistical Computing*. Vienna, Austria: R Foundation for Statistical Computing. <https://www.R-project.org/>.
