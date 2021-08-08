version 15
set matastrict on

mata:
// Function called r_0 in the paper, non-coverage probability
real vector r (real vector t, real vector chi) {
    real vector idx, r, res

    idx = selectindex(sqrt(t) :- chi :<= 5) // Select relevant indices
    r = normal(-sqrt(t):-chi) :+ normal(sqrt(t):-chi)
    res = J(1,length(t), 1)
    if (sum(idx) > 0) res[idx] = r[idx]
    
    return(res)
}


// Derivative of r
real vector r1 (real vector t, real scalar chi) {
    real vector idx, res, case1

    idx = selectindex(t :>= 1e-8) 
    case1 = (normalden(sqrt(t):-chi) :- normalden(sqrt(t):+chi)):/(2*sqrt(t))
    res = J(1,length(t),chi:*normalden(chi))
    if (sum(idx) > 0)  res[idx] = case1[idx]
    
    return(res)
}

// Second Derivative of r
real vector r2 (real vector t, real scalar chi) {
    real vector idx, res, case1

    idx = selectindex(t :>= 2e-6) 
    case1 = (normalden(sqrt(t):+chi):*(chi:*sqrt(t) :+ t :+ 1) :+ ///
        normalden(sqrt(t):-chi):*(chi:*sqrt(t) :- t :- 1)):/(4*t:^(3/2))
    res = J(1,length(t),normalden(chi):*chi:*(chi^2-3)/6)
    if (sum(idx) > 0) res[idx] = case1[idx]
    
    return(res)
}

// Third Derivative of r
real vector r3 (real vector t, real scalar chi) {
    real vector idx, res, case1

    idx = selectindex(t :>= 2e-4) 
    case1 = (normalden(sqrt(t):-chi):*(t:^2 :- 2*chi:*t:^(3/2) :+ (2+chi^2) :*t :- 3*chi:*sqrt(t) :+ 3) :- ///
        normalden(sqrt(t):+chi):*(t:^2 :+ 2*chi:*t:^(3/2) :+ (2+chi^2) :*t :+ 3*chi:*sqrt(t) :+ 3)):/(8*t:^(5/2))
    res = J(1,length(t),normalden(chi):*(chi^5 - 10:*chi^3 + 15:*chi)/60)
    if (sum(idx) > 0) res[idx] = case1[idx]
    
    return(res)
}

// Find t0 and inflection point, called t1 in the paper
struct t0ip {
    real scalar t0, ip
}

real scalar f0 (t,chi) return(r(t,chi) - t*r1(t,chi) - r(0,chi))

struct t0ip scalar rt0 (real scalar chi) {

    transmorphic mroot // avoids printing zeros for mm_root
    struct t0ip scalar res
    real scalar t0, ip, tol, up, lo

    // Find point where we touch origin
    // Need root finding function to be installed
    // ssc install moremata, replace
    tol = 1e-12
    if (chi^2 < 3) {
        t0 = ip = 0
    } 
    else {
        if( abs(r2(chi^2 -3/2, chi)) < tol | (chi^2-3) - chi^2 ==0) ip = chi^2 -3/2 
        else mroot = mm_root(ip=., &r2(), chi^2 -3, chi^2, tol, 1000, chi)

        up = 2*chi^2
        lo = ip
        while (f0(up,chi) <0) {
            lo = up
            up = 2*up
        }
        t0 = lo

        if (f0(lo,chi)<0) {
            mroot = mm_root(t0=., &f0(), lo, up, tol, 1000, chi)
        }
        else if (f0(lo,chi)>tol) printf("{red}Warning: Failed to solve for t0 using rt0 at chi = %9.0g. \n", chi)
    }
    
    res.ip = ip
    res.t0 = t0
    return(res)
}

// \rho(m_2, \chi)
real vector rho0 (real vector t, real vector chi){
    struct t0ip scalar rt0_out
    real scalar res, t0
    real vector idx, case1

    rt0_out = rt0(chi)
    t0 = rt0_out.t0

    idx = selectindex(t :>= t0)
    case1 = r(t,chi)
    res = r(t0,chi) :+ (t:- t0):*r1(t0,chi)
    if (sum(idx) > 0) res[idx] = case1[idx]
    return(res)
}

// \delta(x; x_0)
real vector delta (real vector x, real scalar x0, real scalar chi) {
    real vector res, idx, case2
    idx = selectindex(abs(x:-x0) :> 1e-4)
    res = r2(x0,chi)/2*J(1,length(x),1)
    case2 = (r(x,chi) :- r(x0,chi) :- r1(x0,chi):*(x:-x0)) :/ (x:-x0):^2 
    if (sum(idx) > 0) res[idx] = case2[idx]
    return(res)
}

real vector delta1 (real vector x, real scalar x0, real scalar chi) {
    real vector res, idx, case2

    idx = selectindex(abs(x:-x0) :>= 1e-3)
    res = r3(x0,chi)/6*J(1,length(x),1)
    case2 = ((r1(x,chi) :+ r1(x0,chi)) :- 2*(r(x,chi):-r(x0,chi)):/(x:-x0)) :/ (x:-x0):^2
    if (sum(idx) > 0) res[idx] = case2[idx]
    return(res)
}

// maximize delta(x,x0,chi) over x
struct lamx0 {
    real scalar lam, x0
}

/* Code up evaluator for delta */
struct lamx0 scalar lam (real scalar x0, real scalar chi) {
    struct t0ip scalar rt0_out
    struct lamx0 scalar opt0, res
    real vector xs, rt0_vec, der, val, indDerNeg, iMaxVal, w, diffDerPos, indDerPos
    real vector diffDerPosleq0, iDergeq0, ival, x, rr1max, robj_vec, rmin_vec
    real scalar sumDerNeg, idx, rr1val, i, idxmin, ngrid
    real scalar lo, up, phi, x2,x3, eps, flo, fup, c, xlen, fx2, fx3, fmin 

    // Check derivatives at 0, inflection pt, t0 and x0
    // If above ip, then max is below it, and it's at zero if der at zero is negative
    rt0_out = rt0(chi)
    rt0_vec = (rt0_out.t0, rt0_out.ip)
    xs = (min(rt0_vec), max(rt0_vec))
    if (x0 >= xs[1]) xs = (0,xs[1])
    else {
        xs = uniqrows((0\x0\xs'))
        xs=xs'
    }

    der = delta1(xs,x0,chi)
    val = delta(xs,x0,chi)
    opt0.lam = val[1]; opt0.x0=0
    
    // Expect delta has single maximum, so first increasing, then decr up to tol
    // Need to convert the conditionals to scalars
    indDerNeg = der :<= 0; sumDerNeg = sum(indDerNeg)
    indDerPos = der :>= 0
    maxindex(val,1,iMaxVal,w) // Assign iMaxVal for index of max
    
    diffDerPos = indDerPos[2..length(indDerPos)] - indDerPos[1..(length(indDerPos)-1)]
    diffDerPosleq0 = diffDerPos :<= 0
        
    if (sumDerNeg == length(der) & iMaxVal == 1) {
        return(opt0)
    }
    else if (sum(diffDerPosleq0) == length(diffDerPosleq0) & der[length(der)] <= 0) {
        //Function first increasing, then decreasing
        minindex(indDerPos,1,iDergeq0,w) // Assign iDergeq0 for index of min
        idx = max((min(iDergeq0'),2))
        ival = xs[((idx-1)..idx)]
    }
    else if (min(abs(der)) < 1e-6){
        // Determine interval based on value of delta
        ival = (xs[max((iMaxVal-1,1))],xs[min((iMaxVal+1,length(val)))])
    }
    else {
        //_error("There are multiple optima in the function delta(x, x0, chi)")
        // Commented out warning message
        printf("{red}Warning: There are multiple optima in the function delta(x, x0, chi). \n")
        return(opt0)
    } 

    // Golden search algorithm 
    eps = 1e-8
    lo = min(ival); up = max(ival)
    phi = (sqrt(5) + 1)/2 //golden ratio
    while (up-lo >eps) {
        // Step 1: evaluate end points
        flo = -delta(lo,x0,chi); fup = -delta(up,x0,chi)
        // Step 2: calculate interior point
        c = 2-phi
        xlen = up -lo
        x2 = lo + (xlen)*c
        // Step 3: check convergence criterion
        if (xlen > eps) {
            // Step 4: calculate other interior point
            x3 = up - (xlen)*c
            // Step 5: choose points for next iteration
            fx2 = -delta(x2,x0,chi); fx3 = -delta(x3,x0,chi)
            fmin = min((flo,fx2,fx3,fup))
            if (fmin == flo| fmin == fx2) up =x3 
            else lo = x2
        }
    }
    rr1max = x2
    rr1val = -delta(x2,x0,chi)
    
    
    // Finally check optimum at 0 is not higher
    if (-rr1val > opt0.lam) {
        res.lam = -rr1val; res.x0 = rr1max
        return(res)
    }
    else if (-rr1val > opt0.lam - 1e-9) {
        return(opt0)
    }
    else {
        // Commented out warning message
        printf("{red}Warning: Optimum may be wrong for lam(x0=%9.0g, chi=%9.0g) \n", x0,chi)
        return(opt0)
    }

}

struct alphaxp {
    real vector alpha, x, p
}

// Functions required for rho
real scalar lammax (real scalar x0, real scalar tbar, real scalar chi){
    real scalar res
    struct lamx0 scalar lam_out
    if (x0 >= tbar) {
        res = delta(0,x0,chi)
        return(res)
    }
    else {
        lam_out = lam(x0,chi)
        res = max((lam_out.lam,0))
        return(res)
    }
}

real scalar obj (real scalar x0, real scalar chi, real scalar m2, real scalar tbar, real scalar kappa) {
    real scalar res
    res = r(x0,chi) + r1(x0,chi)*(m2-x0) + lammax(x0,tbar,chi)*(kappa*m2^2-2*x0*m2+x0^2)
    return(res)
}


struct alphaxp scalar rho (real scalar m2, real scalar kappa, real scalar chi) {
    // This function is missing the check for the linear program
    real vector r0, rmin_vec, robj_vec, p, xs0, xs0_vec
    real scalar t0, tbar, i, idxmin, rrmin, rrval, ngrid, rrmina, rrvala, rrminb, rrvalb
    real scalar lo, up, phi, x2,x3, eps, flo, fup, c, xlen, fx2, fx3, fmin 
    struct t0ip scalar rt0_out
    struct alphaxp scalar res
    struct lamx0 scalar lam_out, lam_min

    r0 = rho0(m2,chi)
    rt0_out = rt0(chi)
    t0 = rt0_out.t0

    if (kappa ==1) {
        res.alpha=r(m2,chi); res.x=(0,m2); res.p=(0,1)
        return(res)
    }
    else if (max(abs(m2)) >= t0) {
        // Concave part of parameter space
        res.alpha=r0; res.x=(0,m2); res.p=(0,1)
        return(res)
    }
    else if (kappa >= 1e10 | max(abs(m2:*kappa)) :>= t0) {
        // LF under rho: (0, t0) wp (1-m2/t0, m2/t0), E[t^2]=m2*t0
        // So here kappa doesn't bind
        res.alpha=r0; res.x=(0,t0); res.p=(1:-m2:/t0, m2:/t0)
        return(res)
    }
    else {
        // First determine where delta(x, 0) is maximized
        lam_out = lam(0,chi) 
        tbar = lam_out.x0
                
        // Optimization routine
        //rb
        if (tbar >0) {
            // Golden search
            eps = 1e-8
            lo = 0; up = tbar
            phi = (sqrt(5) + 1)/2 //golden ratio
            while (up-lo >eps) {
                // Step 1: evaluate end points
                flo = obj(lo,chi,m2,tbar,kappa); fup = obj(up,chi,m2,tbar,kappa)
                // Step 2: calculate interior point
                c = 2-phi
                xlen = up -lo
                x2 = lo + (xlen)*c
                // Step 3: check convergence criterion
                if (xlen > eps) {
                    // Step 4: calculate other interior point
                    x3 = up - (xlen)*c
                    // Step 5: choose points for next iteration
                    fx2 = obj(x2,chi,m2,tbar,kappa); fx3 = obj(x3,chi,m2,tbar,kappa)
                    fmin = min((flo,fx2,fx3,fup))
                    if (fmin == flo| fmin == fx2) up =x3 
                    else lo = x2
                }
            }
            rrminb = x2
            rrvalb = obj(x2,chi,m2,tbar,kappa)

        }
        else {
            rrminb = 0
            rrvalb = obj(0,chi,m2,tbar,kappa)
        }

        //ra
        // Golden search
        eps = 1e-8
        lo = tbar; up = t0
        phi = (sqrt(5) + 1)/2 //golden ratio
        while (up-lo >eps) {
            // Step 1: evaluate end points
            flo = obj(lo,chi,m2,tbar,kappa); fup = obj(up,chi,m2,tbar,kappa)
            // Step 2: calculate interior point
            c = 2-phi
            xlen = up -lo
            x2 = lo + (xlen)*c
            // Step 3: check convergence criterion
            if (xlen > eps) {
                // Step 4: calculate other interior point
                x3 = up - (xlen)*c
                // Step 5: choose points for next iteration
                fx2 = obj(x2,chi,m2,tbar,kappa); fx3 = obj(x3,chi,m2,tbar,kappa)
                fmin = min((flo,fx2,fx3,fup))
                if (fmin == flo| fmin == fx2) up =x3 
                else lo = x2
            }
        }
        rrmina = x2
        rrvala = obj(x2,chi,m2,tbar,kappa)
        
        // Comparing ra and rb
        if (rrvalb < rrvala) {
            rrmin = rrminb
            rrval = rrvalb
        }
        else {
            rrmin = rrmina
            rrval = rrvala    
        }
        
        // Post optimization calculations
        lam_min = lam(rrmin,chi)
        xs0_vec = (rrmin, lam_min.x0)
        xs0 = (min(xs0_vec),max(xs0_vec))
        p = (m2-xs0[2])/(xs0[1] - xs0[2])

        res.alpha=rrval; res.x=xs0; res.p=(p, 1-p)
        return(res)
        
    }
}

real vector CVb (real vector B, real scalar alpha) {
    real vector r, idx, case2

    idx = selectindex(B :< 10)
    r = B :+ invnormal(1-alpha)
    case2 = sqrt(invnchi2(1, B:^2,1-alpha))
    if (sum(idx) >0) r[idx] = case2[idx]
    return(r)
}

// Compute average coverage critical value under moment constraints
struct cvaxp {
    real vector cv, alpha, x, p
}
real scalar rho0m_alpha (chi,m2,alpha) return(rho0(m2,chi) - alpha)

real scalar rhom_alpha (chi,m2,kappa,alpha) {
    struct alphaxp scalar rho_out
    rho_out = rho(m2, kappa, chi)
    return(rho_out.alpha - alpha)
}

struct cvaxp scalar cva_fn (real scalar m2, real scalar kappa, real scalar alpha) {
    transmorphic mroot
    struct cvaxp scalar res
    real scalar tol, chi_sol, lo, up, cv, check_val
    real vector limits
    struct alphaxp scalar rho_out

    tol = 1e-12
    if (m2 >= 1e12) {
        res.cv = .; res.alpha = alpha; res.x = (0,m2); res.p=(0,1)
        return(res)
    }
    else if (kappa == 1 | m2==0) {
        res.cv = CVb(sqrt(m2),alpha); res.alpha = alpha; res.x = (0,m2); res.p=(0,1)
        return(res)
    }
    else {
        // Omit case of large m2 for now
        lo = CVb(sqrt(m2),alpha) - 0.01
        limits = (lo,.)
        up = sqrt((1+m2)/alpha)
        if (abs(rho0(m2,up)-alpha) < 9e-6) {
            limits[2] = up
        }
        else {
            mroot = mm_root(chi_sol=., &rho0m_alpha(), lo, up, tol, 1000, m2,alpha)
            limits[2] = chi_sol
        }

        // If rejection rate is already close to alpha, keep cv under kappa = Inf
        check_val = rhom_alpha(limits[2],m2,kappa,alpha)
        if (check_val < -1e-5) {
            mroot = mm_root(cv=., &rhom_alpha(), limits[1], limits[2], tol, 1000, m2,kappa,alpha)
        }
        else {
            cv = limits[2]
        }
        rho_out = rho(m2, kappa, cv)
        res.cv = cv; res.alpha = rho_out.alpha; res.x = rho_out.x; res.p=rho_out.p
        return(res)
        
        
    }
}

// Functions that approximate cva
real scalar f2 (real scalar m2){
    real scalar res
    if (m2 <= 0.002) {
        res = 1.96
    }
    else if (m2 <= 0.4){
        res = -78.639 + 81.46*m2^0.1 - 7.544*log(m2) - 0.26*(log(m2))^2
    }
    else if (m2 <= 10) {
        res = -556.141 + 559.434*m2^0.1 - 54.325*log(m2) - 2.059*(log(m2))^2
    }
    else {
        res = -4.062 + 2.746*m2^0.1 + 4.461*m2^0.5 - 0.582*(log(m2))
    }
    return(res)
}

real scalar cvapprox (real scalar m2, real scalar kappa) {
    real scalar res
    if (m2 > 0.05 & m2 <= 4 & kappa >=1 & kappa <=3) {
        res = -142.179 + 146.044*m2^0.1 -1.417*m2^0.5 + 146.917*kappa^0.1 - 2.339*kappa^0.5 ///
            -146.961*m2^0.1*kappa^0.1 + 2.574*m2^0.5*kappa^0.5 + 0.944*log(m2)*log(kappa)
    }
    else if (m2 > 0.05 & m2 <= 1 & kappa >3 & kappa <=10) {
        res = 14.625 - 18.029*m2^0.1 + 4.915*m2^0.5 - 0.133*(log(m2))^2 ///
            + 0.348*m2^0.5*kappa^0.5 + 1.687 *log(kappa)^(-1) - 0.963*(log(kappa))^(-2)
    }
    else if (m2 > 0.05 & m2 <= 1 & kappa >10 & kappa <=47.2){
        res = 25.496 - 31.484*m2^0.1 + 9.236*m2^0.5 - 0.249*(log(m2))^2 - 0.081*(log(kappa))^(-2)
    }
    else if (m2 > 1 & m2 <= 4 & kappa >3 & kappa <=17.4) {
        res = 26.358 + 2.0495*m2^(0.5) - 7.0787*kappa^0.5 + 2.4885*(log(kappa))^2 - 12.481*m2^0.1*kappa^0.1 ///
            +0.679*m2^0.5*kappa^0.5 + 0.417*log(m2)*log(kappa) - 3.729*(log(kappa))^(-0.5)
    }
    else if (m2 > 4 & m2 <= 10 & kappa >=1 & kappa <=3) {
        res = 15.225 + 1.287*m2^0.5 - 1.339*kappa^0.1 - 12.765*m2^(0.1)*kappa^(0.1) /// 
            + 0.9385*m2^0.5*kappa^0.5 + 0.9125*log(m2)*log(kappa)
    }
    else if (m2 > 4 & m2 <= 10 & kappa >3 & kappa <=10) {
        res = 6.7945 + 0.998*(log(m2))^2 + 15.768*kappa^0.1 - 18.959*m2^0.1*kappa^0.1 ///
            + 0.305*m2^0.5*kappa^0.5 + 1.027*log(m2)*log(kappa)
    }
    else if (m2 > 4 & m2 <= 10 & kappa >10 & kappa <=18.5) {
        res = 88.102 + 245.52*m2^0.1 - 54.885*kappa^0.5 + 24.89*(log(kappa))^2 ///
            - 233.27*m2^0.1*kappa^0.1 + 1.386*m2^0.5*kappa^0.5 + 2.0055*log(m2)*log(kappa)
    }
    else {
        res = f2(m2)
    }
    return(res)
}

real scalar cva (real scalar m2, real scalar kappa, real scalar alpha, real scalar approx) {
    // Choose cva function based on whether to use approx
    struct cvaxp scalar res
    real scalar resapp
    if (approx == 1) {
        resapp = cvapprox(m2, kappa)
        return(resapp)
    }
    else {
        res = cva_fn(m2, kappa,alpha)
        return(res.cv)
    } 
}

real scalar ci_length (real scalar w, real scalar kappa, real scalar alpha, real scalar S, real scalar approx) {
        real scalar m2, res
        m2 = (1/w - 1)^2*S
        res = cva(m2, kappa, alpha, approx)*w
        return(res)
}

real vector Fw_opt (real scalar S, real scalar kappa, real scalar alpha, real scalar approx) {
        real scalar maxbias, minbias, rmin, rval, m2, idxmin, I
        real scalar lo, up, phi, x2,x3, eps, flo, fup, c, xlen, fx2, fx3, fmin 
        real vector robj_vec, rmin_vec, res

        maxbias = 100^2
        minbias = 0
        // Golden search
        eps = 1e-8
        lo = 1/(sqrt(maxbias/S)+1)
        up = 1/(sqrt(minbias/S)+1)
        phi = (sqrt(5) + 1)/2 //golden ratio
        while (up-lo >eps) {
            // Step 1: evaluate end points
            flo = ci_length(lo,kappa,alpha,S, approx); fup = ci_length(up,kappa,alpha,S, approx)
            // Step 2: calculate interior point
            c = 2-phi
            xlen = up -lo
            x2 = lo + (xlen)*c
            // Step 3: check convergence criterion
            if (xlen > eps) {
                // Step 4: calculate other interior point
                x3 = up - (xlen)*c
                // Step 5: choose points for next iteration
                fx2 = ci_length(x2,kappa,alpha,S, approx); fx3 = ci_length(x3,kappa,alpha,S, approx)
                fmin = min((flo,fx2,fx3,fup))
                if (fmin == flo| fmin == fx2) up =x3 
                else lo = x2
            }
        }
        rmin = x2
        rval = ci_length(rmin,kappa,alpha,S, approx)
        
        m2 = (1/rmin - 1)^2*S
        res = J(1,3,.)
        res[1] = rmin // w
        res[2] = cva(m2,kappa,alpha,approx)*rmin // length
        res[3] = m2 //m2
        return(res)
}

real vector Fw_eb (real scalar S, real scalar kappa, real scalar alpha, real scalar approx) {
    real vector res
    real scalar w

    res = J(1,3,.)
    w = S/(1+S)
    res[1] = w //w
    res[2] = cva(1/S,kappa,alpha,approx)*w // length
    res[3] = 1/S //m2
    return(res)
}

// Weighted sample variance and mean of truncated normal
real scalar weighted_mean (real colvector A, real colvector B) {
    // Take in column vectors 
    return(sum(A:*B)/sum(B))
}

real scalar wgtV (real colvector Z, real colvector wgt) {
    return(sum(wgt:^2:*(Z:^2:-weighted_mean(Z, wgt):^2))/(sum(wgt):^2:-sum(wgt:^2)))
}

real scalar tmean (real scalar m, real scalar V) {
    return(m + sqrt(V)*normalden(m/sqrt(V))/normal(m/sqrt(V)))
}

real scalar ExtractAlpha (struct alphaxp scalar rho_out) {
    return(rho_out.alpha)
}

end 

// EBCI regression program
capture program drop ebreg

// Check that moremata is installed and install if not
cap which moremata.hlp
if _rc ssc install moremata

program ebreg, eclass  
// Omit optimal shrinkage for length(se) > 30
// using approx option overwrites alpha option
syntax varlist(min=1 numeric) [if] [in], se(varname) [weights(varname) alpha(real 0.05) kappa(real 1e12)  ///
        wopt approx fs_correction(string) reg_options(string)] 

preserve
marksample touse
qui keep if `touse'

local za invnormal(1-`alpha'/2)
local depvar `1'
local xlist : subinstr local varlist "`depvar'" "", word

if ("`weights'" == "")  reg `varlist', `reg_options'
else reg `varlist' [aw=`weights'], `reg_options'


predict mu1, xb
putmata mu1, replace

putmata se = `se', replace

// Set default method
if ("`fs_correction'" == "") local fs_correction "PMT"

mata{

    Y = st_data(.,"`depvar'")
    za = `za'
    if ("`weights'" == "") wgt = J(length(Y),1,1)
    else wgt = st_data(.,"`weights'")
    approx = "`approx'" == "approx"

    Yt = Y
    set = se
    //delta = st_matrix("e(b)")
    
    w2 = (Yt :- mu1):^2 :- set:^2
    w4 = (Yt :- mu1):^4 :- 6:*set:^2:*(Yt:-mu1):^2 :+ 3:*set:^4
    
    // tilde mu2 and tilde mu4 
    tmu = (weighted_mean(w2,wgt), weighted_mean(w4,wgt))   
    pmt_trim = (2:*mean(wgt:^2:*set:^4) :/ (sum(wgt):*mean(wgt:*set:^2)),
        32:*mean(wgt:^2:*set:^8) :/ (sum(wgt):*mean(wgt:*set:^4)))

    if ("`fs_correction'" == "none") {
        mu2 = max((tmu[1],0))
    }
    else if ("`fs_correction'" == "PMT") {
        mu2 = max((tmu[1],pmt_trim[1]))
    }
    else if ("`fs_correction'" == "FPLIB") {
        mu2 = tmean(tmu[1],wgtV(w2,wgt))
    }

    // Send results back to Stata
    st_numscalar("mu2",mu2)

    if (mu2==0) {
        w_eb = J(length(Yt),1,0)
        w_opt = J(length(Yt),1,0)
        ncov_pa = J(length(Yt),1,.)
        len_eb = J(length(Yt),1,.)
        len_op = J(length(Yt),1,.)
        len_pa = J(length(Yt),1,.)
        len_us = za:*se
        th_us = Y
        th_eb = Y :-mu1
        th_op = J(length(Yt),1,.)
        weights = wgt
        residuals = Y -mu1
        EB_df = (w_eb , w_opt , ncov_pa , len_eb , len_op , len_pa , len_us , th_us , th_eb , th_op, se, weights, residuals)
        st_matrix("EB_df", EB_df)
    }

    else {
        // Assign kappa value
        if (`kappa' == 1e12) {
            if ("`fs_correction'" == "none") {
                kappa = max((tmu[2]/mu2^2,1))
            }
            else if ("`fs_correction'" == "PMT") {
                kappa = max((tmu[2]/mu2^2,1+pmt_trim[2]/mu2^2))
            }
            else if ("`fs_correction'" == "FPLIB") {
                kappa = 1+tmean(tmu[2]-tmu[1]^2,wgtV(w4-2*mu2*w2, wgt))/mu2^2
            }
        }
        else kappa = `kappa'

        n = length(Yt)

        // EB shrinkage
        mu2se = mu2:/(se:^2)
        ebmat = J(n,3,.)
        for (i=1; i <= n; i++) {
            ebmat[i,] = Fw_eb(mu2se[i],kappa,`alpha',approx)
        }
        ebw = ebmat[,1]
        ebl = ebmat[,2]
        ebm = ebmat[,3]
        th_eb = mu1 :+ ebw:*(Y :-mu1)
        w_eb = ebw
        
        
        if ("`wopt'" == "wopt") {
            opmat = J(n,3,.)
            for (i=1; i <= n; i++) {
                opmat[i,] = Fw_opt(mu2se[i],kappa,`alpha',approx)
            }
            opw = opmat[,1]
            opl = opmat[,2]
            opm = opmat[,3]

            th_op = mu1:+opw:*(Yt:-mu1)
            w_opt = opw
        }
        else {
            th_op = J(n,1,.)
            w_opt = J(n,1,.)
            opl = J(n,1,.)
        }

        // Calculate ncov
        ncov_mat = J(n,1,.)
        se2mu = se:^2:/mu2
        chi_vec = za*sqrt(1:+se2mu)
        for (i=1; i <= n; i++) {
            rho_out = rho(se2mu[i],kappa,chi_vec[i])
            ncov_mat[i] = ExtractAlpha(rho_out)
        }
        ncov_pa = ncov_mat
                
        len_eb = ebl:*se
        len_op = opl:*se
        len_pa = sqrt(ebw):*za:*se
        len_us = za:*se
        th_us = Y
        weights = wgt
        residuals = Y -mu1
        
        EB_df = (w_eb , w_opt , ncov_pa , len_eb , len_op , len_pa , len_us , th_us , th_eb , th_op, se, weights, residuals)

        st_numscalar("tmu1",tmu[1])
        st_numscalar("tmu2",tmu[2])
        st_numscalar("kappa",kappa)
        //st_matrix("delta",delta)
        st_matrix("EB_df",EB_df)
    
    }

    // export individual columns of EB output
    st_matrix("w_eb",w_eb)
    st_matrix("w_opt",w_opt)
    st_matrix("ncov_pa",ncov_pa)
    st_matrix("len_eb",len_eb)
    st_matrix("len_op",len_op)
    st_matrix("len_pa",len_pa)
    st_matrix("len_us",len_us)
    st_matrix("th_us",th_us)
    st_matrix("th_eb",th_eb)
    st_matrix("th_op",th_op)
    st_matrix("se",se)
    st_matrix("weights",weights)
    st_matrix("residuals",residuals)
}
restore

// Display results on stata
if (mu2 == 0) {
    di "{red}Warning: mu2 estimate is 0"
    matrix colnames EB_df = "w_eb" "w_opt" "ncov_pa" "len_eb" "len_op" "len_pa" "len_us" "th_us" "th_eb" "th_op" "se" "weights" "residuals"
    di "sqrt_mu2: " 0
    di "kappa = " `kappa'
    di "delta = "
    ereturn matrix EB_df = EB_df
    ereturn scalar mu2 = 0
}
else { 
//getmata (ebw ebl ebm) = ebmat, replace
matrix colnames EB_df ="w_eb" "w_opt" "ncov_pa" "len_eb" "len_op" "len_pa" "len_us" "th_us" "th_eb" "th_op" "se" "weights" "residuals"
di "mu2: estimate = " mu2 ", uncorrected estimate = " tmu1
di "kappa: estimate = " kappa ", uncorrected estimate = " tmu2/tmu1^2
ereturn matrix EB_df = EB_df
ereturn scalar mu2 = mu2
ereturn scalar mu2_uncorrected = tmu1
ereturn scalar kappa = kappa
ereturn scalar kappa_uncorrected = tmu2/tmu1^2
}

ereturn matrix w_eb = w_eb
ereturn matrix w_opt = w_opt
ereturn matrix ncov_pa = ncov_pa
ereturn matrix len_eb = len_eb
ereturn matrix len_op = len_op
ereturn matrix len_pa = len_pa
ereturn matrix len_us = len_us
ereturn matrix th_us = th_us
ereturn matrix th_eb = th_eb
ereturn matrix th_op = th_op
ereturn matrix se = se
ereturn matrix weights = weights
ereturn matrix residuals = residuals

end
