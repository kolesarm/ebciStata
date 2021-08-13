# ebciStata

This Stata package implements robust empirical Bayes confidence intervals from
[Armstrong, Kolesár, and Plagborg-Møller
(2020)](https://arxiv.org/abs/2004.03448). See [ebci_matlab](https://github.com/mikkelpm/ebci_matlab) for a
Matlab version of this package, and [ebci](https://github.com/kolesarm/ebci) for an R version.

## Installation

`ebciStata` is not currently available from SSC. To install directly from this repository, you can copy and run the following lines in Stata:
```stata
// Remove program if it existed previously
cap ado uninstall ebciStata
// Install most up-to-date version
net install ebciStata, from("https://raw.githubusercontent.com/kolesarm/ebciStata/")
```

## Example
```stata
use data/cz, replace
// Use precision weights
gen wgt = 1/se25^2
// Replicate emirical example in Armstrong, Kolesár, and Plagborg-Møller (2020)
ebreg theta25 stayer25, se(se25) wopt weights(wgt)
matrix U = J(rowsof(e(w_opt)),1,1)
matrix means = (U'*e(w_eb), U'*e(w_opt), U'*e(ncov_pa))/e(N)
matrix list means
```
