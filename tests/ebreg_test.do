* Do file to test out function
* Save EB_df to folder to be called by R

clear all
cd ..
use data/cz, replace
gen wgt = 1/se25^2
save tests/output/cz1, replace

ebreg theta25 stayer25 if state == "NY", se(se25) wopt weights(wgt) fs_correction(none)
matrix EB_df = e(EB_df)
svmat EB_df
keep EB_df*
save tests/output/eb_out1, replace

use tests/output/cz1, replace
ebreg theta25 stayer25 if state == "NY", se(se25) wopt weights(wgt) alpha(0.1)
matrix list e(w_opt)
matrix EB_df = e(EB_df)
svmat EB_df
keep EB_df*
save tests/output/eb_out2, replace

use tests/output/cz1, replace
ebreg theta25 stayer25 if state == "FL", se(se25) wopt weights(wgt) kappa(10) alpha(0.1)
matrix EB_df = e(EB_df)
svmat EB_df
keep EB_df*
save tests/output/eb_out3, replace

cd tests/
