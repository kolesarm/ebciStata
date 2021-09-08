**********************
** Empirical Example
**********************

capture log close
log using ebreg_example.log, replace

clear all

** Load Chetty & Hendren (2018) data
* Focus on NY
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

* Compute robust EBCIs for NY
ebreg theta25 $xvar if state == "NY", se(se25) wopt weights(wgt) alpha($alpha) genvar(ebci)
sum ebci_*

* Display various objects obtained
di "lower and upper endpoints of EBCIs for NY"
list ebci_cil_eb ebci_ciu_eb if state == "NY"

di "Shrinkage factors are w_eb and half length len_eb"
list ebci_w_eb ebci_len_eb if state == "NY"

di "Average length of unshrunk CIs relative to robust EBCIs"
qui sum ebci_len_us
local len_us `r(mean)'
qui sum ebci_len_eb
di "Unshrunk CI's are " `len_us'/`r(mean)' " times of robust EBCIs"


********************************************
** Parametric EBCIs (Morris, 1983)
********************************************

* Parametric EBCI results are stored as pa prefix
* Note that the parametric EBCIs are centered at the same shrinkage estimates as the robust EBCIs

di "Average length of parametric EBCIs relative to robust EBCIs"
qui sum ebci_len_pa
local len_pa `r(mean)'
qui sum ebci_len_eb
di "Parametric CI's are " `len_pa'/`r(mean)' " times of robust EBCIs"

di "Average worst-case non-coverage probability of parametric EBCIs"
qui sum ebci_ncov_pa
di "Average worst-case non-coverage probability is " `r(mean)'

* The parametric EBCIs are shorter but less robust: They could potentially
* have Empirical Bayes non-coverage probabilities well above alpha=10%

// Replicate empirical example in Armstrong, Kolesár, and Plagborg-Møller (2020)
/*
ebreg theta25 stayer25, se(se25) wopt weights(wgt) alpha(0.1)
matrix U = J(rowsof(e(w_opt)),1,1)
matrix means = U'*(e(w_eb), e(w_opt), e(ncov_pa), e(len_eb), e(len_op), e(len_pa), e(len_us))/e(N)
matrix list means
*/


log close
