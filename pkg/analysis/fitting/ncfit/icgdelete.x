include	<gset.h>
include	<mach.h>
include	<pkg/gtools.h>
include	"icfit.h"

define	MSIZE		3.		# Mark size

# ICG_DELETE -- Delete data point nearest the cursor.
# The nearest point to the cursor in NDC coordinates is determined.

procedure icg_delete (ic, gp, gt, nl, x, y, wts, userwts, npts, wx, wy)

pointer	ic					# ICFIT pointer
pointer	gp					# GIO pointer
pointer	gt					# GTOOLS pointer
pointer	nl					# NLFIT pointer
real	x[npts], y[npts]			# Data points
real	wts[npts], userwts[npts]		# Weight arrays
int	npts					# Number of points
real	wx, wy					# Position to be nearest
#--

int	gt_geti()
pointer	sp, xout, yout

begin
	call smark (sp)
	call salloc (xout, npts, TY_REAL)
	call salloc (yout, npts, TY_REAL)

	call icg_axes (ic, gt, nl, 1, x, y, Memr[xout], npts)
	call icg_axes (ic, gt, nl, 2, x, y, Memr[yout], npts)

	if (gt_geti (gt, GTTRANSPOSE) == NO)
	    call icg_d1r (ic, gp, gt, Memr[xout], Memr[yout], wts, userwts,
		npts, wx, wy)
	else
	    call icg_d1r (ic, gp, gt, Memr[yout], Memr[xout], wts, userwts,
		npts, wy, wx)

	call sfree (sp)
end


# ICG_D1 -- Do the actual delete.

procedure icg_d1r (ic, gp, gt, x, y, wts, userwts, npts, wx, wy)

pointer	ic					# ICFIT pointer
pointer	gp					# GIO pointer
pointer	gt					# GTOOLS pointer
real	x[npts], y[npts]			# Data points
real	wts[npts], userwts[npts]		# Weight arrays
int	npts					# Number of points
real	wx, wy					# Position to be nearest

int	i, j
real	x0, y0, r2, r2min, aux

begin
	# Transform world cursor coordinates to NDC.

	call gctran (gp, wx, wy, wx, wy, 1, 0)

	# Search for nearest point to a point with non-zero weight.

	r2min = MAX_REAL
	do i = 1, npts {
	    if (wts[i] == 0.)
		next

	    call gctran (gp, real (x[i]), real (y[i]), x0, y0, 1, 0)
		
	    r2 = (x0 - wx) ** 2 + (y0 - wy) ** 2

	    # Whith mag scale, uses only x dimension. 
	    if (IC_YAXIS(ic) == IC_MAG) {
	        if (abs(x0 - wx) < r2min) {
		    r2min = abs (x0 - wx)
		    j = i
	        }
	    } else {
	        if (r2 < r2min) {
		    r2min = r2
		    j = i
	        }
	    }
	}

	# Mark the deleted point with a cross and set the weight to zero.

	if (j != 0) {
	    aux = real (y[j])
	    if (IC_YAXIS(ic) == IC_MAG) {
	        if (aux > IC_MINY(ic))
	            aux = IC_MAG0(ic) - 2.5 * log10 (aux)
	        else
	            aux = IC_MAG0(ic) - 2.5 * log10 (IC_MINY(ic))
	    }
	    call gscur (gp, real (x[j]), aux)
	    call gmark (gp, real (x[j]), aux, GM_CROSS, MSIZE, MSIZE)
	    wts[j] = 0.
	    IC_NEWWTS(ic) = YES
	}
end
