include <pkg/gtools.h>
include	<error.h>
include	<mach.h>
include	<ctype.h>
include <fset.h>	# FIO
include	<gset.h>	# GIO

include "../ncfit/names.h"		# Ensures use of correct icfit 


# NLFIT1 -- Given the x, y and error data determine the 1-d fitting function.
# If the interactive flag is set then set the fitting parameters interactively.
#
# Modified 2/22/96 (I Busko) to implement chi-squared minimization.

procedure nl_fit1 (x, y, err, npix, nl, ic, gt, filename, interactive, 
		   graphics, out)

real	x[ARB]				# input x data
real	y[ARB]				# input y data
real	err[ARB]			# input err data
int	npix				# number of data points
pointer	ic				# ICFIT (modified) pointer
pointer	gt				# GTOOLS pointer
pointer	nl				# NLFIT pointer
char	filename[ARB]			# original file name
bool	interactive			# Interactive?
char	graphics[SZ_FNAME]		# Graphics device
pointer	out				# Table pointer

#--
int compare()
extern	compare
int	i, nl_units()
pointer	gp, sp, ttx, tx, ty, twts, wts, indpt, gopen()
real	dummy[1]

int	icg_fit()

common	/comm/ttx

errchk	gopen, icg_fit, ic_fit

begin
	# Make sure INDEF-valued coefficients are replaced
	# by something sensible.
	call nl_startv (nl, x, dummy, y)

	# Allocate memory for curve fitting.
	call smark (sp)
	call salloc (wts, npix, TY_REAL)
	call salloc (tx , npix, TY_REAL)
	ttx	= tx
	call salloc (ty ,   npix, TY_REAL)
	call salloc (twts , npix, TY_REAL)
	call salloc (indpt, npix, TY_INT)

	# Prepare weight vector.
	call prep_data (x, y, err, npix, nl, wts)


	# Sort the data for use before fitting.
	do i = 1, npix
	   Memi[indpt+i-1] = i
	call amovr( x, Memr[tx], npix)
	call tqsort ( Memi[indpt], npix, compare)
	do i = 1, npix {
	    Memr[tx+i-1] = x[ Memi[indpt+i-1] ]
	    Memr[ty+i-1] = y[ Memi[indpt+i-1] ]
	    Memr[twts+i-1] = Memr[wts+ Memi[indpt+i-1] - 1]
	}

	call ic_putr (ic, "xmin", Memr[tx])
	call ic_putr (ic, "xmax", Memr[tx+npix-1])

	# Check x axis and x ref. units. 
	if (nl_units (nl, x, npix) == ERR)
	    call error (0, "Impossible to plot data and function.")

	# If the interactive flag is set then use icg_fit to set the
	# fitting parameters.  
	# The weights are reset since the user may delete points

	if (interactive) {
	    gp = gopen (graphics, NEW_FILE, STDGRAPH)

	    i = icg_fit (ic, gp, "gcur", gt, nl, 
			 Memr[tx], Memr[ty], Memr[twts], npix)
	    call gclose (gp)
	} else {
	    call ic_fit (ic, nl, Memr[tx], Memr[ty], Memr[twts], npix,
			 YES, YES, YES, YES)
	    i = OK
	}

	if (i == OK)
	    call nl_wtdb( filename, nl, out, true )
	else
	    call eprintf ("No data written to output table.\n")

	call sfree (sp)
end
                                                          
                                                              

                                                       
