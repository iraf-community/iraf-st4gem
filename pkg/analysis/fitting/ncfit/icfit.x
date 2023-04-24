include	"icfit.h"

# IC_FIT -- Fit a function.  This is the main fitting task.

procedure ic_fit (ic, nl, x, y, wts, npts, newx, newy, newwts, newfunction)

pointer	ic				# ICFIT pointer
pointer	nl				# NLFIT pointer
real	x[npts]				# Ordinates
real	y[npts]				# Data to be fit
real	wts[npts]			# Weights
int	npts				# Number of points
int	newx				# New x points?
int	newy				# New y points?
int	newwts				# New weights?
int	newfunction			# New function?
#--

pointer	sp, dummy

errchk	ic_dosetup, ic_reject, nl_fit, rg_wtbinr

begin
	call smark (sp)
	call salloc (dummy, npts, TY_REAL)
	call amovkr (0., Memr[dummy], npts)

	# Setup the new parameters.
	call ic_dosetup (ic, nl, x, wts, npts, newx, newwts, newfunction)

	# If necessary, clear previously rejected data points.
	if ((IC_LOW(ic) == 0.) && (IC_HIGH(ic) == 0.)) 
	    call nl_clear (nl)
	if (IC_NITERATE(ic) == 0)
	    call nl_clear (nl)

	# If not sampling use the data array directly.
	if (npts == IC_NFIT(ic)) {
	        call nl_puti (nl, "npts", npts)
	        call nl_fit (nl, x, Memr[dummy], y, wts)

	# If sampling first form the sample y values.
	} else {
	    if ((newx == YES) || (newy == YES) || (newwts == YES))
	        call rg_wtbinr (IC_RG(ic), IC_NAVERAGE(ic), y, wts, npts,
		    Memr[IC_YFIT(ic)], Memr[IC_WTSFIT(ic)], IC_NFIT(ic))
	    # Update number of data points and fit.
	    call nl_puti (nl, "npts", IC_NFIT(ic))
	    call nl_fit (nl, Memr[IC_XFIT(ic)], Memr[dummy], Memr[IC_YFIT(ic)],
			 Memr[IC_WTSFIT(ic)])
	}

	# Do pixel rejection if desired.
	if ((IC_LOW(ic) > 0.) || (IC_HIGH(ic) > 0.)) {
	    if (npts == IC_NFIT(ic)) {
	        call ic_reject (nl, x, y, wts, Memi[IC_REJPTS(ic)],
		    IC_NFIT(ic), IC_LOW(ic), IC_HIGH(ic), IC_NITERATE(ic),
		    IC_GROW(ic), IC_NREJECT(ic))
	    } else {
	        call ic_reject (nl, Memr[IC_XFIT(ic)], Memr[IC_YFIT(ic)],
		    Memr[IC_WTSFIT(ic)], Memi[IC_REJPTS(ic)], IC_NFIT(ic),
		    IC_LOW(ic), IC_HIGH(ic), IC_NITERATE(ic), IC_GROW(ic),
		    IC_NREJECT(ic))
	    }
	} else
	    IC_NREJECT(ic) = 0

	call sfree (sp)
end
