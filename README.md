--- 
title: "Replicating Ledoit and Wolf (2008, 2011) Functions" 
date: '`r format(Sys.time(), "%d, %B, %Y")`'
author: "Paulo Ferreira Naibert" 
output: 
  html_document:
    toc: false
    number_sections: true
keep_tex: true
bibliography: refs.bib
--- 


# Introduction

On @lw2008-sr and @lw2011-var, the authors provide a method to test the hypothesis of difference in the Sharpe Ratio and the Variance of two investment strategies.
The papers can be found in the following links:

[SR Paper](https://www.econ.uzh.ch/dam/jcr:ffffffff-935a-b0d6-0000-00007214c2bc/jef_2008pdf.pdf), 
[SR WP](http://www.econ.uzh.ch/static/wp_iew/iewwp320.pdf), 
[Variance Paper](https://www.econ.uzh.ch/dam/jcr:520edf26-2322-4708-8dde-51d61141914a/Ledoit_et_al-2011-Wilmott_Robust_Performance.pdf), 
[Variance WP](http://www.econ.uzh.ch/static/wp_iew/iewwp516.pdf).

In this post we are going to replicate @lw2008-sr and @lw2011-var.
The codes can be found in 
[Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html), 
in the 
[publications section](https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html).
Our goal is to verify if I replicated the functions on Michael Wolf's website.
On what follows is a test of the functions' outputs.


# System Information

## Software Information
```{r}
sessionInfo()
```

# DISCLAIMER

I own NONE of the rights to the codes.
I simply edited the already existing codes in
[Michael Wolf's Webpage](https://www.econ.uzh.ch/en/people/faculty/wolf.html).
I also offer NO support for the functions.
They are presented AS IS.

# REFERENCES