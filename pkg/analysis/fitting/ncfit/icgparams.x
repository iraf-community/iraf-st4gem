include	<pkg/gtools.h>
include	"icfit.h"

# ICG_PARAMS -- Set parameter string.

procedure icg_params (ic, nl, x, y, wts, npts, gt)

pointer	ic		# ICFIT pointer
pointer	nl		# NLFIT pointer
real	x[ARB]		# Ordinates
real	y[ARB]		# Abscissas
real	wts[ARB]	# Weights
int	npts		# Number of data points
pointer	gt		# GTOOLS pointer
#--

int	i, n, deleted
real	rms, chisq
pointer	sp, fit, wts1, str, params, dummy

real	nl_statr()

begin
	call smark (sp)

	n = IC_NFIT(ic)
	deleted = 0
	rms   = INDEFR
	chisq = INDEFR

	if (n == npts) {
	    # Allocate memory for the fit.

	    call salloc (fit, n, TY_REAL)
	    call salloc (wts1, n, TY_REAL)
	    call salloc (dummy, n, TY_REAL)
	    call amovkr (0., Memr[dummy], n)

	    # Eliminate rejected points and count deleted points.

	    call amovr (wts, Memr[wts1], n)
	    if (IC_NREJECT(ic) > 0) {
		do i = 1, npts {
		    if (Memi[IC_REJPTS(ic)+i-1] == YES)
			Memr[wts1+i-1] = 0.
		}
	    }
	    deleted = 0
	    do i = 1, n {
		if (wts[i] == 0.)
		    deleted = deleted + 1
	    }

	    # Set the fit and compute the RMS error.

	    if (IC_FITERROR(ic) == NO) {
	        call nl_vector (nl, x, Memr[dummy], Memr[fit], n)
	        rms   = nl_statr (nl, "rms")
	        chisq = nl_statr (nl, "chisq")
	    } else {
		rms   = INDEFR
		chisq = INDEFR
	    }
	} else if (n > 0) {
	    # Allocate memory for the fit.

	    call salloc (fit, n, TY_REAL)
	    call salloc (wts1, n, TY_REAL)
	    call salloc (dummy, n, TY_REAL)	# added by PEH 21-May-1992
            call amovkr (0., Memr[dummy], n)

	    # Eliminate rejected points and count deleted points.

	    call amovr (Memr[IC_WTSFIT(ic)], Memr[wts1], n)
	    if (IC_NREJECT(ic) > 0) {
	        do i = 1, npts {
		    if (Memi[IC_REJPTS(ic)+i-1] == YES)
			Memr[wts1+i-1] = 0.
		}
	    }
	    deleted = 0
	    do i = 1, n {
		if (wts[i] == 0.)
		    deleted = deleted + 1
	    }

	    # Set the fit and compute the rms error.

	    if (IC_FITERROR(ic) == NO) {
	        call nl_vector (nl, Memr[IC_XFIT(ic)], Memr[dummy], 
				Memr[fit], n)
		rms   = nl_statr (nl, "rms")
		chisq = nl_statr (nl, "chisq")
	    } else {
		rms   = INDEFR
		chisq = INDEFR
	    }
	}

	# Print the parameters and errors.

	call salloc (str, SZ_LINE, TY_CHAR)
	call salloc (params, 2*SZ_LINE, TY_CHAR)

	call sprintf (Memc[str], SZ_LINE,
	    "func=%s, low_rej=%r, high_rej=%r, niterate=%d, grow=%r")
	    call ic_gstr (ic, "function", Memc[params], 2*SZ_LINE)
	    call pargstr (Memc[params])
	    call pargr (IC_LOW(ic))
	    call pargr (IC_HIGH(ic))
	    call pargi (IC_NITERATE(ic))
	    call pargr (IC_GROW(ic))

	call sprintf (Memc[params], 2*SZ_LINE,
	    "%s\ntotal=%d, sample=%d, rejected=%d, deleted=%d, CHISQ=%7.4g")
	    call pargstr (Memc[str])
	    call pargi (npts)
	    call pargi (n)
	    call pargi (IC_NREJECT(ic))
	    call pargi (deleted)
#	    call pargr (rms)
	    call pargr (chisq)

	call gt_sets (gt, GTPARAMS, Memc[params])

	# Free allocated memory.

	call sfree (sp)
end
