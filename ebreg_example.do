** Do file to test out function
capture log close
log using ebreg_example.log, replace

clear all
use data/cz, replace
gen wgt = 1/se25^2

** Show how to replicate first two columns of Table 1


log close
