# ebciStata

This Stata package implements robust empirical Bayes confidence intervals from
[Armstrong, Kolesár, and Plagborg-Møller
(2021)](https://arxiv.org/abs/2004.03448). See
[ebci_matlab](https://github.com/mikkelpm/ebci_matlab) for a Matlab version of
this package, and [ebci](https://github.com/kolesarm/ebci) for an R version.

This software package is based upon work supported by the National Science
Foundation under grant numbers SES-2049765 (Armstrong), SES-22049356 (Kolesár),
and SES-1851665 (Plagborg-Møller), and by work supported by the Alfred P. Sloan
Research Fellowship (Kolesár).

## Installation

`ebciStata` is not currently available from SSC. To install directly from this
repository, you can copy and run the following lines in Stata:
```stata
// Remove program if it existed previously
capture ado uninstall ebreg
// Install most up-to-date version
net install ebreg, from("https://raw.githubusercontent.com/kolesarm/ebciStata/master")
```

## Example

Estimates and robust EBCIs for neighborhood effects, as in the empirical
application in [Armstrong, Kolesár, and Plagborg-Møller
(2021)](https://arxiv.org/abs/2004.03448). Shrink fixed-effect estimates of the
neighborhood effects, for children with parents at the 25th percentile of the
income distribution (`theta25`) toward average outcome for permanent residents
(stayers) at the 25th percentile of the income distribution. Use precision
weights proportional to the inverse of the squared standard error of the
fixed-effect estimates.

```stata
use data/cz, clear
// Use precision weights
generate wgt = 1/se25^2
ebreg theta25 stayer25, se(se25) weights(wgt) alpha(0.1) genvar(ebci)
/* Shrinkage estimates, neighborhood effects, and confidence intervals
   are stored in variables with prefix specified in the genvar option.
   Summarize them.
*/
summarize ebci_*

// List EB point estimate and robust EBCI for NY commuting zones
list czname ebci_th_eb ebci_cil_eb ebci_ciu_eb if state=="NY"
```
