Replication of Ledoit and Wolf (2008, 2011) Functions
================
Paulo Ferreira Naibert

Original: 2019 October 13.

Current: 16, outubro, 2019

Overview
========

This Repository contains functions to replicate Ledoit and Wolf (2008) and Ledoit and Wolf (2011).

Michael Wolf provides `R` codes for those two papers on [his webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html), in the [publications section](https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html).

What is found here is merely an edit of those codes.

PAPERS and GOALS of the code
============================

On Ledoit and Wolf (2008) and Ledoit and Wolf (2011), the authors provide methods to test the hypothesis of difference in the Sharpe Ratio and the Variance of two investment strategies. The papers can be found in the following links:

[SR Paper](https://www.econ.uzh.ch/dam/jcr:ffffffff-935a-b0d6-0000-00007214c2bc/jef_2008pdf.pdf), [SR WP](http://www.econ.uzh.ch/static/wp_iew/iewwp320.pdf), [Variance Paper](https://www.econ.uzh.ch/dam/jcr:520edf26-2322-4708-8dde-51d61141914a/Ledoit_et_al-2011-Wilmott_Robust_Performance.pdf), [Variance WP](http://www.econ.uzh.ch/static/wp_iew/iewwp516.pdf).

Functions
=========

I copied the original functions on `Sharpe.RData` and `Var.RData` into the files: `"SR-ORIGINAL.R"` and `"Var-ORIGINAL.R"` with some minor modifications. Then I copied the functions to the file `"LW-SR-Var.R"` and modified them more.

I didn't have enought time to modify the `bloc.size.calibrate()` functions.

Data
====

For the data, I simpy saved the returns data of `Sharpe.RData` as a `.RDS` file. Loading `rets.RDS` should result in a list object on the `R` Software (R Core Team 2019) .

TODO
====

-   Document Functions (.tex?)
-   Add technincal notes
-   Add Mean test
-   Add E-V utility test (West, Edison, and Cho 1993, Kirby and Ostdiek (2012))

DISCLAIMER
==========

I own NONE of the rights to the codes. I simply edited the already existing codes in [Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html). I also offer NO support for the functions. They are presented "AS IS".

REFERENCES
==========

Kirby, Chris, and Barbara Ostdiek. 2012. “It’s All in the Timing: Simple Active Portfolio Strategies That Outperform Naïve Diversification.” *Journal of Financial and Quantitative Analysis* 47 (2): 437–67.

Ledoit, Oliver, and Michael Wolf. 2008. “Robust performance hypothesis testing with the Sharpe ratio.” *Journal of Empirical Finance* 15 (5): 850–59.

———. 2011. “Robust performance hypothesis testing with the Variance.” *Wilmott Magazine*, no. 55 (September): 86–89.

R Core Team. 2019. *R: A Language and Environment for Statistical Computing*. Vienna, Austria: R Foundation for Statistical Computing. <https://www.R-project.org/>.

West, Kenneth D., Hali J. Edison, and Dongchul Cho. 1993. “A utility-based comparison of some models of exchange rate volatility.” *Journal of International Economics* 35 (1-2): 23–45.
