include	"../nlfit/nlfit.h"

# PREP_DATA -- Prepare the x,y,err data to be fitted.
#
# This routine looks for problems in the x,y,err vectors as INDEF-valued
# data or errors, zeroed errors and such, and cleans up to the standards
# needed by the nlfit/curfit packages. This means setting the weight vector
# properly, as well as setting eventual x,y weird values to something safe. 
# Routine also builds the weights to allow chi-square computation.
#
# Created 2/21/96 (I Busko).

procedure prep_data (x, y, err, npix, nl, wts)

real	x[ARB]		# io: x data
real	y[ARB]		# io: y data
real	err[ARB]	# i:  err data
int	npix		# i:  number of data points
pointer	nl		# i:  NLFIT pointer
pointer	wts		# o:  weight vector.

int	i
real	xsafe, ysafe
real	sigma, epadu, rnoise, rn2

int	nl_stati()
real	nl_statr()

begin
	sigma  = nl_statr (nl, "sigma")
	epadu  = nl_statr (nl, "epadu")
	rnoise = nl_statr (nl, "readnoise")
	rn2 = rnoise * rnoise

	do i = 1, npix {

	    # If any INDEF value, or zeroed error bar, set weight to zero 
	    # once and for all. Otherwise, error is a valid one, thus build
	    # weight for chi-square computation. 
	    if ( IS_INDEFR (x[i])   || 
	         IS_INDEFR (y[i])   || 
	         IS_INDEFR (err[i]) ||
	         (err[i] == 0.0)    ) {
	        Memr[wts+i-1] = 0.0
	    } else {
		Memr[wts+i-1] = 1.0 / err[i]

	        # If sigma parameter was provided by user, use it  
	        # instead of any error bars present in the data.
	        if (!IS_INDEF(sigma))
	            Memr[wts+i-1] = 1.0 / sigma

	        # If ccd type was set, this takes preference.
	        if (nl_stati(nl,"errtyp") == CCD) {
	            if (y[i] <= 0.) {
	                if (rnoise > 0)
	                    Memr[wts+i-1] = 1.0 / rnoise
	                else
	                    Memr[wts+i-1] = 0.0
	            } else
	                Memr[wts+i-1] = 1.0 / sqrt (y[i]/epadu + rn2)
	        }
	    }
	}

	# Set all the zero weighted values to something safe.
	do i = 1, npix {
	    if ( Memr[wts+i-1] != 0.0 ) {
		xsafe = x[i]
		ysafe = y[i]
	    }
	}
	do i = 1, npix {
	    if ( Memr[wts+i-1] == 0.0) {
		x[i] = xsafe
		y[i] = ysafe
	    }
	}

end
