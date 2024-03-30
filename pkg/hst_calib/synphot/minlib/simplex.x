#* HISTORY *
#* B.Simon	12-Sep-94	Original

# SIMPLEX -- Minimize function by downhill simplex method

procedure simplex (fitfunc, ftol, maxiter, nprint, ndata, nvar, p, y)

extern	fitfunc		# i: function to fit
double	ftol		# i: termination condition
int	maxiter		# i: maximum number of iterations
int	nprint		# i: number of iterations between diagnostic prints
int	ndata		# i: number of data points
int	nvar		# i: number of free variables in fit
double	p[nvar+1,nvar]	# u: simplex containing trial solution
double	y[nvar+1]	# u: function values at simplex points
#--
double	epsil, rtol, xtol, delta, ypr, yprr, alpha, beta, gamma
int	iter, iflag, mpts, ilo, ihi, inhi, i, j
pointer	sp, pr, prr, pbar, resid, errmsg

data	alpha, beta, gamma  / 1.0, 0.5, 2.0 /
string	baditer   "Fit terminated without convergence"

double	dpmpar(), enorm()

begin
	call smark (sp)
	call salloc (pr, nvar, TY_DOUBLE)
	call salloc (prr, nvar, TY_DOUBLE)
	call salloc (pbar, nvar, TY_DOUBLE)
	call salloc (resid, ndata, TY_DOUBLE)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	iter = 0
	mpts = nvar + 1
        epsil = dpmpar (1)
	rtol = max (epsil, ftol)

	repeat {
	    # Find lowest (ilo), highest (ihi) and next highest (inhi)
	    # points in the simplex

	    ilo = 1
	    if (y[1] > y[2]) {
		ihi = 1
		inhi = 2
	    } else {
		ihi = 2
		inhi = 1
	    }

	    do i = 1, mpts {
		if (y[i] < y[ilo])
		    ilo = i

		if (y[i] > y[ihi]) {
		    inhi = ihi
		    ihi = i
		} else if (y[i] > y[inhi]) {
		    if (i != ihi)
			inhi = i
		}
	    }

	    # Check for ending conditions

	    xtol = 0.0
	    do j = 1, nvar {
		delta = 0.5 * (p[ihi,j] + p[ilo,j])
		if (delta < epsil)
		    delta = max (abs(p[ilo,j]), epsil)

		delta = (p[ihi,j] - p[ilo,j]) / delta
		xtol = xtol + delta * delta
	    }

	    if (sqrt(xtol) <= rtol)
		break

	    if (iter == maxiter) {
		call sprintf (Memc[errmsg], SZ_FNAME, "%d")
		call pargi (maxiter)

		call synphotwarn (baditer, Memc[errmsg])
		break
	    }

	    # Compute average of all points except highest, i.e. the center 
	    # of the face of the simplex opposite the highest point. We will
	    # explore along the ray joining these two points to find a new
	    # point for the simplex.

	    iter = iter + 1
	    do j = 1, nvar 
		Memd[pbar+j-1] = 0.0

	    do i = 1, mpts {
		if (i != ihi) {
		    do j = 1, nvar
			Memd[pbar+j-1] = Memd[pbar+j-1] + p[i,j]
		}
	    }

	    # Reflect the simplex around the highest point

	    do j = 1, nvar {
		Memd[pbar+j-1] = Memd[pbar+j-1] / nvar
		Memd[pr+j-1] = (1.0 + alpha) * Memd[pbar+j-1] - 
			       alpha * p[ihi,j]
	    }

	    iflag = 1
            call fitfunc (ndata, nvar, Memd[pr], Memd[resid], iflag)
	    ypr = enorm (ndata, Memd[resid])
	    if (iflag < 0)
		break

	    # If reflection gives a better solution, try an additional
	    # extrapolation in the same direction

	    if (ypr <= y[ilo]) {
		do j = 1, nvar 
		    Memd[prr+j-1] = gamma * Memd[pr+j-1] + 
				    (1.0 - gamma) * Memd[pbar+j-1]

		iflag = 2
		call fitfunc (ndata, nvar, Memd[prr], Memd[resid], iflag)
		yprr = enorm (ndata, Memd[resid])
		if (iflag < 0)
		    break

		if (yprr < y[ilo]) {
		    do j = 1, nvar
			p[ihi,j] = Memd[prr+j-1]
		    y[ihi] = yprr
		} else {
		    do j = 1, nvar
			p[ihi,j] = Memd[pr+j-1]
		    y[ihi] = ypr
		}

	    } else if (ypr >= y[inhi]) {
		# The reflected point is worse than the second highest
		# Check to see if it worse than the highest. If it's
		# better, replace the highest.

		if (ypr < y[ihi]) {
		    do j = 1, nvar
			p[ihi,j] = Memd[pr+j-1]
		    y[ihi] = ypr
		}

		# But look for an intermediate lowest point

		do j = 1, nvar 
		    Memd[prr+j-1] = beta * p[ihi,j] + 
				    (1.0 - beta) * Memd[pbar+j-1]

		iflag = 2
		call fitfunc (ndata, nvar, Memd[prr], Memd[resid], iflag)
		yprr = enorm (ndata, Memd[resid])
		if (iflag < 0)
		    break

		if (yprr < y[ihi]) {
		    do j = 1, nvar
			p[ihi,j] = Memd[prr+j-1]
		    y[ihi] = yprr

		} else {
		    # Can't find a point better than the highest,
		    # so contract simplex around the lowest point

		    do i = 1, mpts {
			if (i != ilo) {
			    do j = 1, nvar {
				Memd[pr+j-1] = 0.5 * (p[i,j] + p[ilo,j])
				p[i,j] = Memd[pr+j-1]
			    }

			    iflag = 2
			    call fitfunc (ndata, nvar, Memd[pr], 
					  Memd[resid], iflag)
			    y[i] = enorm (ndata, Memd[resid])
			}
		    }
		}

	    } else {
		# If the original reflection gives a middling point,
		# replace the highest point

		do j = 1,nvar
		    p[ihi,j] = Memd[pr+j-1]
		y[ihi] = ypr
	    }

            # Diagnostic print of intermediate results

	    if (nprint != 0) {
		iflag = 0
		if (mod (iter-1, nprint) == 0)
		    call fitfunc (ndata, nvar, Memd[pr], Memd[resid], iflag)
	    }
	}

	call sfree (sp)
end
