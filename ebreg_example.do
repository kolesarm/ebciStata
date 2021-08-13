** Replicate first two columns of Table 1
capture log close
log using ebreg_example.log, replace

clear all
use data/cz, replace
// Use precision weights
gen wgt = 1/se25^2

// Replicate emirical example in Armstrong, Kolesár, and Plagborg-Møller (2020)
ebreg theta25 stayer25, se(se25) wopt weights(wgt) alpha(0.1)
matrix U = J(rowsof(e(w_opt)),1,1)
matrix means = U'*(e(w_eb), e(w_opt), e(ncov_pa), e(len_eb), e(len_op), e(len_pa), e(len_us))/e(N)
matrix list means

log close
