**********************
** Empirical Example
**********************

capture log close
log using ebreg_example.log, replace

clear all

** Load Chetty & Hendren (2018) data
* Display results for NY
use data/cz, replace

* Preliminary estimates of causal effect Y_i: fixed effects estimates of
* neighborhood effect, for children with parents at the 25th percentile 
* of the income distribution
list theta25 se25 if state == "NY"

* Shrink toward average outcome for permanent residents (stayers) at the
* 25th percentile of the income distribution
global xvar stayer25

* Precision weights for regression and moment estimates
gen wgt = 1/se25^2

* Significance level
global alpha 0.1

********************************************
** Robust EBCIs, baseline shrinkage
********************************************

* Compute robust EBCIs for whole population, display results for NY
ebreg theta25 $xvar, se(se25) wopt weights(wgt) alpha($alpha) genvar(ebci)

* Replicate empirical example in Armstrong, Kolesár, and Plagborg-Møller (2020)
sum ebci_*

* Display various objects obtained
di "lower and upper endpoints of EBCIs for NY"
list ebci_cil_eb ebci_ciu_eb if state == "NY"

di "Shrinkage factors are w_eb and half length len_eb"
list ebci_w_eb ebci_len_eb if state == "NY"

di "Average length of unshrunk CIs relative to robust EBCIs"
qui sum ebci_len_us if state == "NY"
local len_us `r(mean)'
qui sum ebci_len_eb if state == "NY"
di "Unshrunk CI's are " `len_us'/`r(mean)' " times of robust EBCIs"


********************************************
** Parametric EBCIs (Morris, 1983)
********************************************

* Parametric EBCI results are stored as pa prefix
* Note that the parametric EBCIs are centered at the same shrinkage estimates as the robust EBCIs

di "Average length of parametric EBCIs relative to robust EBCIs"
qui sum ebci_len_pa if state == "NY"
local len_pa `r(mean)'
qui sum ebci_len_eb if state == "NY"
di "Parametric CI's are " `len_pa'/`r(mean)' " times of robust EBCIs"

di "Average worst-case non-coverage probability of parametric EBCIs"
qui sum ebci_ncov_pa if state == "NY"
di "Average worst-case non-coverage probability is " `r(mean)'

* The parametric EBCIs are shorter but less robust: They could potentially
* have Empirical Bayes non-coverage probabilities well above alpha=10%


log close
