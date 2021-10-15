{smcl}
{* *! version 1.0.0 15oct2021}{...}
{title:Title}
{viewerjumpto "Syntax" "ebreg##syntax"}{...}
{viewerjumpto "Description" "ebreg##description"}{...}
{viewerjumpto "Examples" "ebreg##examples"}{...}
{viewerjumpto "Stored Results" "ebreg##results"}{...}
{viewerjumpto "References" "ebreg##references"}{...}
{pstd}
{hi:ebreg} {hline 2} Robust Empirical Bayes Confidence Intervals


{marker syntax}{title:Syntax}

{p 8 16 2}
{cmd:ebreg} {it:depvar [indepvars] [if] [in]} {cmd:,} {cmd:se(se_var)} {it:[options]}
{p_end}

{synoptset 17}{...}
{synopthdr}
{synoptline}
{synopt :{opt weights(wgt_var)}}An optional numeric vector of weights to be used in the fitting process in computing {it:delta}, {it:mu_2} and {it:kappa}.{p_end}
{synopt :{opt alpha(#)}}Determines confidence level, {cmd:1-alpha}. Default is {cmd:alpha=0.05}.{p_end}
{synopt :{opt kappa(#)}}If specified, use this value for the kurtosis {it:kappa} of {it:theta-X*delta} (such as setting it to infinity), instead of computing it.{p_end}
{synopt :{opt wopt}}Compute length-optimal robust robust empirical Bayes confidence intervals (EBCIs).
These are robust EBCIs centered at estimates with the shrinkage factor chosen to minimize the length of the resulting EBCI.{p_end}
{synopt :{opt fs_correction()}}Finite-sample correction method used to compute mu_2 and kappa.
These corrections ensure that we do not shrink the preliminary estimates (as specified by {cmd:depvar}) all the way to the regression line.
If "PMT", use posterior mean truncation, if "FPLIB" use limited information Bayesian approach with a flat prior,
and if "none", truncate the estimate of mu_2 at 0 and the estimate of kappa at 1. Default is "PMT".{p_end}
{synopt :{opt reg_options(#)}}Passes options into Stata's default regression.{p_end}
{synopt :{opt approx}}Use a polynomial approximation to the critical value (for faster computation), rather than computing it exactly. Assumes {cmd:alpha=0.05}. User-specified option for {cmd:alpha} is ignored. {p_end}
{synopt :{opt genvar(ebci)}}Store {cmd: ebreg} output as variables using prefix ebci_*. The stored output vectors are described below.{p_end}
{synoptline}
{p2colreset}{...}

{it:se()} must be specified. They are standard errors {it:sigma} associated with the preliminary estimates as specified by {cmd:depvar}.


{marker description}{...}
{title:Description}

{pstd}
This Stata package implements robust empirical Bayes confidence intervals from Armstrong, Kolesár, and Plagborg-Møller (2021) by shrinking preliminary estimates toward a regression line.


{marker examples}{...}
{title:Examples}

{phang}{cmd:. use data/cz, clear}{p_end}
{phang}{cmd:. gen wgt = 1/se25^2}{p_end}
{phang}{cmd:. ebreg theta25 stayer25, se(se25) weights(wgt) alpha(0.1) genvar(ebci)}{p_end}
{phang}{cmd:. summarize ebci_*}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ebreg} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(mu2)}}Estimated second moment of {it:theta-X*delta}, {it:mu_2}, computed using the finite-sample correction as specified by {cmd:fs_correction}.{p_end}
{synopt:{cmd:e(mu2_uncorrected)}}Estimated second moment of {it:theta-X*delta}, {it:mu_2} without a finite-sample correction.{p_end}
{synopt:{cmd:e(kappa)}}Estimated kurtosis {it:kappa} of {it:theta-X*delta}, computed using the finite-sample correction as specified by {cmd:fs_correction}.{p_end}
{synopt:{cmd:e(kappa_uncorrected)}}Estimated kurtosis {it:kappa} of {it:theta-X*delta}  without a finite-sample correction.{p_end}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(th_op)}}Vector of estimates based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(th_eb)}}Vector of empirical Bayes estimates.{p_end}
{synopt:{cmd:e(th_us)}}Vector of unshrunk estimates (as specified by depvar).{p_end}
{synopt:{cmd:e(ciu_us)}}Vector of upper endpoints of unshrunk CI.{p_end}
{synopt:{cmd:e(cil_us)}}Vector of lower endpoints of unshrunk CI.{p_end}
{synopt:{cmd:e(ciu_pa)}}Vector of upper endpoints of parametric EBCIs.{p_end}
{synopt:{cmd:e(cil_pa)}}Vector of lower endpoints of parametric EBCIs.{p_end}
{synopt:{cmd:e(ciu_op)}}Vector of upper endpoints of robust EBCIs based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(cil_op)}}Vector of lower endpoints of robust EBCIs based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(ciu_eb)}}Vector of upper endpoints of robust EBCIs.{p_end}
{synopt:{cmd:e(cil_eb)}}Vector of lower endpoints of robust EBCIs.{p_end}
{synopt:{cmd:e(len_us)}}Vector of half-lengths of unshrunk CIs.{p_end}
{synopt:{cmd:e(len_pa)}}Vector of half-lengths of parametric EBCIs.{p_end}
{synopt:{cmd:e(len_op)}}Vector of half-lengths of robust EBCIs based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(len_eb)}}Vector of half-lengths of robust EBCIs.{p_end}
{synopt:{cmd:e(ncov_pa)}}Vector of maximal non-coverages of parametric EBCIs.{p_end}
{synopt:{cmd:e(w_opt)}}Vector of length-optimal shrinkage factors.{p_end}
{synopt:{cmd:e(w_eb)}}Vector of empirical Bayes shrinkage factors, {it:mu_2/(mu_2+sigma^2_i)}.{p_end}
{synopt:{cmd:e(EB_df)}}Matrix of all ebreg output.{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{phang}
Armstrong, Timothy B., Kolesár, Michal, and Plagborg-Møller, Mikkel (2021): Robust Empirical Bayes Confidence Intervals, {browse "https://arxiv.org/abs/2004.03448"}

{phang}You can find a MATLAB version of this command at {browse "https://github.com/mikkelpm/ebci_matlab"}, and an R version at {browse "https://github.com/kolesarm/ebci"}
{p_end}
