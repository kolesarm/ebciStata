{smcl}
{* *! version 1.0.0 Oct 2021}{...}
{title:Title}

{pstd}
{hi:ebreg} {hline 2} Robust Empirical Bayes Confidence Intervals


{title:Syntax}

{p 8 16 2}
{cmd:ebreg} {it:depvar [indepvars] [if] [in]} {cmd:,} {cmd:se(se_var)} {it:[options]}
{p_end}

{synoptset 17}{...}
{synopthdr}
{synoptline}
{synopt :{opt weights(wgt_var)}}An optional numeric vector of weights to be used in the fitting process in computing delta, mu_2 and kappa.{p_end}
{synopt :{opt alpha(#)}}Determines confidence level, 1-alpha. Default is alpha=0.05.{p_end}
{synopt :{opt kappa(#)}}If specified, use this value for the kurtosis kappa of theta-X*delta (such as setting it to infinity), instead of computing it.{p_end}
{synopt :{opt wopt}}Compute length-optimal robust robust empirical Bayes confidence intervals (EBCIs). These are robust EBCIs centered at estimates with the shrinkage factor chosen to minimize the length of the resulting EBCI.{p_end}
{synopt :{opt fs_correction()}}Finite-sample correction method used to compute mu_2 and kappa. These corrections ensure that we do not shrink the preliminary estimates Y all the way to zero. If "PMT", use posterior mean truncation, if "FPLIB" use limited information Bayesian approach with a flat prior, and if "none", truncate the estimates at 0 for mu_2 and 1 for kappa. Default is "PMT".{p_end}
{synopt :{opt reg_options(#)}}Passes options into Stata's default regression.{p_end}
{synopt :{opt approx}}Use a polynomial approximation to the critical value (for faster computation), rather than computing cva() exactly. Assumes alpha=0.05. User-specified option for "alpha" is ignored. {p_end}
{synopt :{opt genvar(ebci)}}Store {cmd: ebreg} output as variables using prefix ebci_*. Output stored are w_eb w_opt ncov_pa len_eb len_op len_pa len_us cil_eb ciu_eb
    cil_op ciu_op cil_pa ciu_pa cil_us ciu_us th_us th_eb th_op described below.{p_end}
{synoptline}
{p2colreset}{...}

{it:se()} must be specified. They are standard errors {it:sigma} associated with the preliminary estimates {it:Y}.


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
{synopt:{cmd:e(mu2)}}Estimated second moment of theta-X*delta, mu_2. Estimate after the finite-sample correction as specified by fs_correction.{p_end}
{synopt:{cmd:e(mu2_uncorrected)}}Estimated second moment of theta-X*delta, mu_2 without fs_correction.{p_end}
{synopt:{cmd:e(kappa)}}Estimated kurtosis kappa of theta-X*delta with fs_correction.{p_end}
{synopt:{cmd:e(kappa_uncorrected)}}Estimated kurtosis kappa of theta-X*delta without fs_correction.{p_end}
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
{synopt:{cmd:e(w_eb)}}Vector of empirical Bayes shrinkage factors, mu_2/(mu_2+sigma^2_i).{p_end}
{synopt:{cmd:e(EB_df)}}Matrix of all ebreg output.{p_end}
{p2colreset}{...}
