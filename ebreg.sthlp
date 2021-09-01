{smcl}
{* *! version 1.0.0 Feb 2021}{...}
{title:Title}

{pstd}
{hi:ebreg} {hline 2} Empirical Bayes Regression.


{title:Syntax}

{p 8 16 2}
{cmd:ebreg} {it:depvar [indepvars] [if] [in]} {cmd:,} {cmd:se(se_var)} {it:[options]}
{p_end}

{synoptset 17}{...}
{synopthdr}
{synoptline}
{synopt :{opt weights(wgt_var)}}An optional vector of weights to be used in the fitting process in computing delta, mu_2 and kappa. Should be NULL or a numeric vector.{p_end}
{synopt :{opt alpha(#)}}Determines confidence level, 1-alpha.{p_end}
{synopt :{opt kappa(#)}}If non-NULL, use pre-specified value for the kurtosis kappa of theta-X*delta (such as Inf), instead of computing it.{p_end}
{synopt :{opt wopt}}Compute length-optimal robust EBCIs. These are robust EBCIs centered at length-optimal shrinkage factors.{p_end}
{synopt :{opt fs_correction()}}Finite-sample correction method used to compute mu_2 and kappa. These corrections ensure that we do not shrink the preliminary estimates Y all the way to zero.
If "PMT", use posterior mean truncation, if "FPLIB" use limited information Bayesian approach with a flat prior, and if "none",
truncate the estimates at 0 for mu_2 and 1 for kappa. Default is "PMT".{p_end}
{synopt :{opt reg_options(#)}}Passes options into Stata's default regression.{p_end}
{synopt :{opt approx}}Use a polynomial approximation to the critical value (for faster computation), rather than computing cva() exactly. 
Assumes alpha=0.05. User-specified option for "alpha" is ignored. {p_end}
{synopt :{opt genvar(ebci)}}Store EBCI output as variables using prefix ebci_*. Output stored are w_eb w_opt ncov_pa len_eb len_op len_pa len_us cil_eb ciu_eb 
    cil_op ciu_op cil_pa ciu_pa cil_us ciu_us th_us th_eb th_op described below.{p_end}
{synoptline}
{p2colreset}{...}

{it:se()} must be specified. They are standard errors {it:sigma} associated with the preliminary estimates {it:Y}.


{marker description}{...}
{title:Description}

{pstd}
This Stata package implements robust empirical Bayes confidence intervals from Armstrong, Kolesár, and Plagborg-Møller (2020).


{marker examples}{...}
{title:Examples}

{phang}{cmd:. use cz, replace}{p_end}
{phang}{cmd:. gen wgt = 1/se25^2}{p_end}
{phang}{cmd:. ebreg theta25 stayer25 if state == "NY", se(se25) wopt weights(wgt) fs_correction(none) }{p_end}
{phang}{cmd:. matrix list e(w_opt)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ebreg} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(mu2)}}Scalar of estimated second moment of theta-X*delta, mu_2^(1/2). Estimate after the finite-sample correction as specified by fs_correction.{p_end}
{synopt:{cmd:e(mu2_uncorrected)}}Scalar of estimated second moment of theta-X*delta, mu_2^(1/2) without fs_correction.{p_end}
{synopt:{cmd:e(kappa)}}Scalar of Estimated kurtosis kappa of theta-X*delta with fs_correction.{p_end}
{synopt:{cmd:e(kappa_uncorrected)}}Scalar of Estimated kurtosis kappa of theta-X*delta without fs_correction.{p_end}
{synopt:{cmd:e(EB_df)}}Matrix of ebreg output{p_end}
{synopt:{cmd:e(residuals)}}Matrix of regression residuals{p_end}
{synopt:{cmd:e(weights)}}Matrix of weights used{p_end}
{synopt:{cmd:e(se)}}Matrix of Standard error sigma, as supplied by the argument se.{p_end}
{synopt:{cmd:e(th_op)}}Matrix of Estimate based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(th_eb)}}Matrix of EB estimate.{p_end}
{synopt:{cmd:e(th_us)}}Matrix of Unshrunk estimate Y{p_end}
{synopt:{cmd:e(len_us)}}Matrix of Half-length of unshrunk CIs{p_end}
{synopt:{cmd:e(len_pa)}}Matrix of Half-length of parametric EBCIs{p_end}
{synopt:{cmd:e(len_op)}}Matrix of Half-length of robust EBCIs based on length-optimal shrinkage{p_end}
{synopt:{cmd:e(cil_eb)}}Matrix of CI lower bound of robust EBCIs based on EB shrinkage{p_end}
{synopt:{cmd:e(ciu_eb)}}Matrix of CI upper bound of robust EBCIs based on EB shrinkage{p_end}
{synopt:{cmd:e(cil_op)}}Matrix of CI lower bound of robust EBCIs based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(ciu_op)}}Matrix of CI upper bound of robust EBCIs based on length-optimal shrinkage.{p_end}
{synopt:{cmd:e(cil_pa)}}Matrix of CI lower bound of robust EBCIs based on parametric EBCIs.{p_end}
{synopt:{cmd:e(ciu_pa)}}Matrix of CI upper bound of robust EBCIs based on parametric EBCIs.{p_end}
{synopt:{cmd:e(cil_us)}}Matrix of unshrunk CI lower bound.{p_end}
{synopt:{cmd:e(ciu_us)}}Matrix of unshrunk CI upper bound.{p_end}
{synopt:{cmd:e(ncov_pa)}}Matrix of Maximal non-coverage of parametric EBCIs{p_end}
{synopt:{cmd:e(w_opt)}}Matrix of Optimal shrinkage factors{p_end}
{synopt:{cmd:e(w_eb)}}Matrix of EB shrinkage factors, mu_2/(mu_2+sigma^2_i){p_end}
{p2colreset}{...}
