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
use cz, replace
gen wgt = 1/se25^2
ebreg theta25 stayer25 if state == "NY", se(se25) wopt weights(wgt) fs_correction(none)
matrix list e(w_opt)
```
