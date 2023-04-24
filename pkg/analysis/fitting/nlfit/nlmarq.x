include "nlfit.h"

# Macros for accessing array elements.

define	I1	Memi[$1+($2-1)]			# 1-D int
define	B1	Memb[$1+($2-1)]			# 1-D bool
define	A1	Memr[$1+($2-1)]			# 1-D real
define	A2	Memr[$1+($3-1)*ncoeff+($2-1)]	# 2-D real
define	D1	Memd[$1+($2-1)]			# 1-D double


#  NL_MARQ --  Levenberg-Marquardt minimization method.
#
#  This module is based on the FORTRAN version of Numerical Recipes. 
#  Several modifications were introduced, such as:
#
#  - stopping criterion takes lambda value into account.
#  - dynamic memory allocation.
#  - layering on top of nlfit data structures and function calls.
#  - specific function computation routine instead of external.
#  - single-size covariance and curvature matrices.
#  - error handling.
#
#  On input, array NL_SPARAMS(nl) contains an initial guess, and on 
#  return the fitted coefficients. On return, array NL_PERRORS(nl)
#  stores the coefficient errors derived from the covariance matrix.
#  On return, NL_CHISQ(nl) contains the chi-squared, and NL_RMS(nl) 
#  contains the r.m.s residual of the fit. 
#
#  On input, array NL_PFLAGS(nl) must be set in order to flag (true) 
#  coefficients which are allowed to vary during the fit. 
#
#  In case a singular matrix is found, ERR is returned and a warning 
#  message is issued, otherwise OK is returned.
#

int procedure nl_marq (nl)

pointer	nl		# i:  Curve descriptor.
#--
pointer	sp
pointer	alpha, covar
int	ncoeff, iter, maxit, i
real	lambda, chisq, ochisq, rms, hold

int	nl_stati(), nl_mqmin()

errchk	nl_min

begin
	ncoeff = nl_stati (nl, "npar")
	maxit  = nl_stati (nl, "maxit")

	# Alloc work space.
	call smark (sp)
	call salloc (alpha, ncoeff * ncoeff, TY_REAL)
	call salloc (covar, ncoeff * ncoeff, TY_REAL)

	# Initialize.
	lambda = -1.
	if (nl_mqmin (nl, covar, alpha, chisq, rms, lambda) == ERR) {
	    call sfree (sp)
	    return (ERR)
	}

	# Iterate until no more improvement in chi-squared or until
	# step size becomes negligible.
	iter = 0
	while (iter < maxit) {
	    ochisq = chisq
	    if (nl_mqmin (nl, covar, alpha, chisq, rms, lambda) == ERR) {
	        call sfree (sp)
	        return (ERR)
	    }
	    if (((chisq < ochisq)              && 
                 (abs(chisq - ochisq) < 0.01)) ||
                (lambda > 1.E30))
	        break
	    iter = iter + 1
	}

	# Finish.
	lambda = 0.
	if (nl_mqmin (nl, covar, alpha, chisq, rms, lambda) == ERR) {
	    call sfree (sp)
	    return (ERR)
	}

	call nl_putr (nl, "chisq", chisq)
	call nl_putr (nl, "rms", rms)

	# Errors are taken directly from the covariance matrix.
	do i = 1, ncoeff {
	    hold = A2(covar,i,i)
	    if (hold > 0.0)
	        D1(NL_PERRORS(nl),i) = sqrt (double(hold))
	    else
	        D1(NL_PERRORS(nl),i) = 0.0D0
	}

	# If reached maximum number of iterations, issue msg.
	if (iter == maxit) {
	    call eprintf("Stop after %d iterations. Chi-sq = %g\n")
	        call pargi (iter)
	        call pargr (chisq)
	    call flush (STDERR)
	}
	call sfree (sp)
	return (OK)
end



#  NL_MQMIN --  Perform one Levenberg-Marquardt iteration.
#
#  "lambda" is the fudge factor for the scale length constant. If
#  set to any negative value, routine initializes the trial coefficient
#  vector and the covariance and curvature matrices. After convergency,
#  this routine must be called with lambda = 0 so the final covariance
#  and chisq can be returned. 
#
#  Matrices "alpha" and "covar" must be allocated by the caller.

int procedure nl_mqmin (nl, covar, alpha, chisq, rms, lambda)

pointer	nl		# i:  Curve descriptor.
pointer	alpha		# io: Curvature matrix.
pointer	covar		# io: Covariance matrix.
real	lambda		# io: Scale length fudge factor.
real	chisq		# o:  Chi-squared.
real	rms		# o:  Rms residual.
#--
pointer	sp, atry, beta, da
int	j, k, l, ncoeff, mfit
real	ochisq

int	nl_stati(), nl_gj()

errchk	nl_mcof, nl_gj, nl_covs

begin
	ncoeff = nl_stati (nl, "npar")

	# Alloc work arrays.
	call smark (sp)
	call salloc (atry, ncoeff, TY_REAL)
	call salloc (beta, ncoeff, TY_REAL)
	call salloc (da,   ncoeff, TY_REAL)

	# Initialization.
	if (lambda < 0.) {
	    mfit = 0
	    do j = 1, ncoeff {
	        if (B1(NL_PFLAGS(nl),j)) 
	            mfit = mfit + 1
	    }
	    lambda = 0.001
	    call nl_mcof (nl, Memr[NL_SPARAMS(nl)], ncoeff, mfit, alpha,
                          Memr[beta], chisq, rms)
	    ochisq = chisq
	    do j = 1, ncoeff 
	        A1(atry,j) = A1(NL_SPARAMS(nl),j)
	}

	# Change fitting matrix by scaling diagonal.
	do j = 1, mfit {
	    do k = 1, mfit
	        A2(covar,j,k) = A2(alpha,j,k)
	    A2(covar,j,j) = A2(alpha,j,j) * (1. + lambda)
	    A1(da,j) = A1(beta,j)
	}

	# Solve.
	if (nl_gj (covar, mfit, ncoeff, Memr[da]) == ERR) {
	    call sfree (sp)
	    return (ERR)
	}

	# Compute covariance matrix at convergency.
	if (lambda == 0.) {
	    call nl_cvsrt (nl, covar, ncoeff, mfit)
	    call sfree (sp)
	    return (OK)
	}

	# See if trial succeeded.
	j = 0
	do l = 1, ncoeff {
	    if (B1(NL_PFLAGS(nl),l)) {
	        j = j + 1
	        A1(atry,l) = A1(NL_SPARAMS(nl),l) + A1(da,j)
	    }
	}
	call nl_mcof (nl, Memr[atry], ncoeff, mfit, covar, Memr[da], 
                     chisq, rms)

	# Yes. Accept new solution.
	if (chisq < ochisq) {
	    lambda = 0.1 * lambda
	    ochisq = chisq
	    do j = 1, mfit {
	        do k = 1, mfit
	            A2(alpha,j,k) = A2(covar,j,k)
	        A1(beta,j) = A1(da,j)
	    }
	    do l = 1 ,ncoeff
	        A1(NL_SPARAMS(nl),l) = A1(atry,l)

	# No. Decrease step and return.
	} else {
	    lambda = 10. * lambda
	    chisq = ochisq
	}

	call sfree (sp)
	return (OK)
end



#  NL_MCOF --  Compute curvature matrix and chi-squared gradient vector.
#
#  This routine also updates the chi-squared and rms values for the current
#  coefficients. It calls routine nl_fdev which computes both the function
#  and its derivates respect each coefficient for a given value of the
#  independent variable(s).

procedure nl_mcof (nl, coeffic, ncoeff, mfit, alpha, beta, chisq, rms)

pointer	nl		# i:  Curve descriptor.
real	coeffic[ARB]	# i:  Coefficient array.
int	ncoeff		# i:  Total number of coefficients.
int	mfit		# i:  Number of fitted coefficients.
pointer	alpha		# io: Curvature matrix.
real	beta[ARB]	# io: Gradient vector.
real	chisq		# o:  Chi-squared.
real	rms		# o:  Rms residual.
#--
pointer	sp, deriv
int	i, j, k, l, m, npts, ndata
real	dy, dy2, sig2i, wt, ymod

int	nl_stati()

errchk	nl_fdev

begin
	# Alloc vector for derivatives.
	call smark (sp)
	call salloc (deriv, ncoeff, TY_REAL)

	# Initialize.
	do j = 1, mfit {
	    do  k = 1, j
	        A2(alpha,j,k) = 0.
	}
	call amovkr (0.0, beta, mfit)

	# Sum over data arrays.
	chisq = 0.
	rms = 0.
	npts = 0
	ndata = nl_stati (nl, "npts")
	do i = 1, ndata {

	    # Skip rejected/flagged data points.
	    if ( !(B1(NL_REJFLAG(nl),i))) {

	        # Compute function and derivatives.
	        npts = npts + 1
	        call nl_fdev (nl, A1(NL_SX(nl),i), A1(NL_SY(nl),i),
                              coeffic, ymod, Memr[deriv])

	        # Compute residuals.
	        sig2i = A1(NL_SW(nl),i) * A1(NL_SW(nl),i)
	        dy    = A1(NL_SZ(nl),i) - ymod
	        dy2   = dy * dy
	        j = 0

	        # Update curvature and gradient arrays.
	        do l = 1, ncoeff {

	            # Skip frozen coefficients.
	            if (B1(NL_PFLAGS(nl),l)) {

	                j = j + 1
	                wt = A1(deriv,l) * sig2i
	                k = 0

	                # Update symmetric curvature matrix.
	                do m = 1, l {
	                    if (B1(NL_PFLAGS(nl),m)) {
	                        k = k + 1
	                        A2(alpha,j,k) = A2(alpha,j,k) + wt * A1(deriv,m)
	                    }
	                }

	                # Update gradient.
	                beta[j] = beta[j] + dy * wt
	            }
	        }

	        # Update stats.
	        chisq = chisq + dy2 * sig2i
	        rms   = rms   + dy2
	    }
	}

	# Compute chidq and rms.
	if (npts > mfit) {
	    chisq = chisq / (npts - mfit)
	    if (rms > 0.)
	        rms   = sqrt (rms / (npts - mfit))
	}

	# Build symmetric matrix.
	do j = 2, mfit {
	    do k = 1, j-1
	        A2(alpha,k,j) = A2(alpha,j,k)
	}

	call sfree (sp)
end



#  NL_CVSRT --  Rearrange covariance matrix into the full "ncoeff"
#               space. Zeros are left in places corresponding to
#               frozen coefficients.

procedure nl_cvsrt (nl, covar, ncoeff, mfit)

pointer	nl		# i:  Curve descriptor.
pointer	covar		# io: Covariance matrix.
int	ncoeff		# i:  Total number of coefficients.
int	mfit		# i:  Number of fitted coefficients.
#--
int	i,j,k
real	swap

begin
	do i  = mfit + 1, ncoeff {
	    do j = 1, i {
	        A2(covar,i,j) = 0.
	        A2(covar,j,i) = 0.
	    }
	}
	k = mfit
	do j = ncoeff, 1, -1 {
	    if (B1(NL_PFLAGS(nl),j)) {
	        do i = 1, ncoeff {
	            swap          = A2(covar,i,k)
	            A2(covar,i,k) = A2(covar,i,j)
	            A2(covar,i,j) = swap
	        }
	        do i = 1, ncoeff {
	            swap          = A2(covar,k,i)
	            A2(covar,k,i) = A2(covar,j,i)
	            A2(covar,j,i) = swap
	        }
	        k = k - 1
	    }
	}
end



#  NL_GJ --  Solve linear system by Gauss-Jordan elimination. If singular
#            matrix, return ERR.

int procedure nl_gj (a, n, ncoeff, b)

pointer	a		# io: System matrix; on output, its inverse.
int	n		# i:  Its size.
int	ncoeff		# i:  Total size allocated by caller (used by macros)
real	b[ARB]		# io: Right-side vector, on output, solution.
#--
pointer	sp, indxc, indxr, ipiv
int	i, icol, irow, j, k, l, ll
real	big, dum, pivinv

begin
	# Alloc and initailize arrays for pivot bookeeping.
	call smark (sp)
	call salloc (indxc, ncoeff, TY_INT)
	call salloc (indxr, ncoeff, TY_INT)
	call salloc (ipiv,  ncoeff, TY_INT)
	call amovki (0, Memi[ipiv], 50)

	# Main loop over columns.
	do i = 1, n {
	    big = 0.

	    # Search for pivot.
	    do j = 1, n {
	        if (I1(ipiv,j) != 1) {
	            do k = 1, n {
	                if (I1(ipiv,k) == 0) {
	                    if (abs (A2(a,j,k)) >= big) {
	                        big = abs (A2(a,j,k))
	                        irow = j
	                        icol = k
	                    }
	                } else if (I1(ipiv,k) > 1) {
	                    call eprintf ("Singular matrix in nl_gj\n")
	                    call sfree (sp)
	                    return (ERR)
	                }
	            }
	        }
	    }
	    I1(ipiv,icol) = I1(ipiv,icol) + 1

	    # Found pivot; interchange rows until pivot is in diagonal.
	    # Use index arrays to avoid physically moving data; just
            # re-label the columns.
	    if (irow != icol) {
	        do l = 1, n {
	            dum          = A2(a,irow,l)
	            A2(a,irow,l) = A2(a,icol,l)
	            A2(a,icol,l) = dum
	        }
	        dum     = b[irow]
	        b[irow] = b[icol]
	        b[icol] = dum
	    }

	    # Divide pivot row by pivot element.
	    I1(indxr,i) = irow
	    I1(indxc,i) = icol
	    if (A2(a,icol,icol) == 0.0D0) {
	        call eprintf ("2 Singular matrix in nl_gj\n")
	        call sfree (sp)
	        return (ERR)
	    }
	    pivinv = 1.0 / A2(a,icol,icol)
	    A2(a,icol,icol) = 1.0D0
	    do l = 1, n
	        A2(a,icol,l) = A2(a,icol,l) * pivinv
	    b[icol] = b[icol] * pivinv

	    # Now reduce the rows.
	    do ll = 1, n {
	        if(ll != icol) {
	            dum = A2(a,ll,icol)
	            A2(a,ll,icol) = 0.
	            do l = 1, n
	                A2(a,ll,l) = A2(a,ll,l) - A2(a,icol,l) * dum
	            b[ll] = b[ll] - b[icol] * dum
	        }
	    }
	}

	# Unscramble the inverse matrix.
	do l = n, 1, -1 {
	    if (I1(indxr,l) != I1(indxc,l)) {
	        do k = 1, n {
	            dum                 = A2(a,k,I1(indxr,l))
	            A2(a,k,I1(indxr,l)) = A2(a,k,I1(indxc,l))
	            A2(a,k,I1(indxc,l)) = dum
	        }
	    }
	}
	call sfree (sp)
	return (OK)
end

