include	"icfit.h"

# IC_DOSETUP -- Setup the fit.  
# This is called at the start of each call to ic_fit to update the fitting 
# parameters if necessary.

procedure ic_dosetup (ic, nl, x, wts, npts, newx, newwts, newfunction)

pointer	ic				# ICFIT pointer
pointer	nl				# NLFIT pointer
real	x[npts]				# Ordinates of data
real	wts[npts]			# Weights
int	npts				# Number of points in data
int	newx				# New x points?
int	newwts				# New weights?
int	newfunction			# New function?
#--

pointer	rg_xrangesr()

errchk	malloc, rg_xrangesr

begin
	# Set sample points.
	if ((newx == YES) || (newwts == YES)) {
	    if (npts == 0)
		call error (0, "No data points for fit")

	    call mfree (IC_XFIT(ic), TY_REAL)
	    call mfree (IC_YFIT(ic), TY_REAL)
	    call malloc (IC_XFIT(ic), npts, TY_REAL)

	    call mfree (IC_WTSFIT(ic), TY_REAL)
	    call malloc (IC_WTSFIT(ic), npts, TY_REAL)

	    call mfree (IC_REJPTS(ic), TY_INT)
	    call malloc (IC_REJPTS(ic), npts, TY_INT)
	    call amovki (NO, Memi[IC_REJPTS(ic)], npts)
	    IC_NREJECT(ic) = 0
	    call nl_clear (nl)

	    # Set sample points.
	    call rg_free (IC_RG(ic))
	    IC_RG(ic) = rg_xrangesr (Memc[IC_SAMPLE(ic)], x, npts)
	    call rg_order (IC_RG(ic))
	    call rg_merge (IC_RG(ic))
	    call rg_wtbinr (IC_RG(ic), max (1, abs (IC_NAVERAGE(ic))), x, wts,
		npts, Memr[IC_XFIT(ic)], Memr[IC_WTSFIT(ic)], IC_NFIT(ic))

	    if (IC_NFIT(ic) == 0)
		call error (0, "No sample points for fit")

	    if (IC_NFIT(ic) == npts) {
	        call rg_free (IC_RG(ic))
	        call mfree (IC_XFIT(ic), TY_REAL)
	        call mfree (IC_WTSFIT(ic), TY_REAL)
	        IC_YFIT(ic) = NULL
		IC_WTSFIT(ic) = NULL
	    } else
	        call malloc (IC_YFIT(ic), IC_NFIT(ic), TY_REAL)
	}
end
