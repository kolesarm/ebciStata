* Do file to test out function
capture log close
log using ebci_example_log, replace

clear all
use cz, replace
gen wgt = 1/se25^2

ebreg theta25 stayer25 if state == "NY", se(se25) wopt weights(wgt) fs_correction(none) 

ebreg theta25 stayer25 if state == "NY", se(se25) wopt weights(wgt) alpha(0.1)
matrix list e(w_opt)

ebreg theta25 stayer25 if state == "FL", se(se25) wopt weights(wgt) kappa(10) alpha(0.1)

ereturn list
display e(kappa)
matrix list e(EB_df)
matrix list e(w_opt)

log close
