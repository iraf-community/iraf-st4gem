include	<pkg/gtools.h>
include	"icfit.h"

# ICG_AXES -- Set axes data.
# The applications program may set additional axes types.

procedure icg_axes (ic, gt, nl, axis, x, y, z, npts)

pointer	ic				# ICFIT pointer
pointer	gt				# GTOOLS pointer
pointer	nl				# NLFIT pointer
int	axis				# Output axis
real	x[npts]				# Independent variable
real	y[npts]				# Dependent variable
real	z[npts]				# Output values
int	npts				# Number of points
#--

int	i, axistype, gtlabel[2], gtunits[2]
real 	a, b, xmin, xmax
pointer	label, units

int	nl_stati()
real	nl_zeval()
errchk	adivr()
real    icg_dvzr() # Error action to take on zero division
extern	icg_dvzr

data	gtlabel/GTXLABEL, GTYLABEL/
data	gtunits/GTXUNITS, GTYUNITS/

begin
	i = nl_stati (nl, "units")
	switch (i) {
	case ANGSTROM:
	    call strcpy ("Wavelength", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("Angstrom", Memc[IC_UNITS(ic,1)], SZ_LINE)
	case CM:
	    call strcpy ("Wavelength", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("cm", Memc[IC_UNITS(ic,1)], SZ_LINE)
	case METER:
	    call strcpy ("Wavelength", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("meter", Memc[IC_UNITS(ic,1)], SZ_LINE)
	case HZ:
	    call strcpy ("Frequency", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("Hz", Memc[IC_UNITS(ic,1)], SZ_LINE)
	case KEV:
	    call strcpy ("Energy", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("KeV", Memc[IC_UNITS(ic,1)], SZ_LINE)
	case LINEAR:
	    call strcpy ("Distance", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("", Memc[IC_UNITS(ic,1)], SZ_LINE)
	case FOURTH_ROOT:
	    call strcpy ("Distance**1/4", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("", Memc[IC_UNITS(ic,1)], SZ_LINE)
	default:
	    call strcpy ("X", Memc[IC_LABELS(ic,1)], SZ_LINE)
	    call strcpy ("", Memc[IC_UNITS(ic,1)], SZ_LINE)
	}

	axistype = IC_AXES(ic, IC_GKEY(ic), axis)
	switch (axistype) {
	case 'x':	# Independent variable
	    call gt_sets (gt, gtlabel[axis], Memc[IC_LABELS(ic,1)])
	    call gt_sets (gt, gtunits[axis], Memc[IC_UNITS(ic,1)])
	    call amovr (x, z, npts)
	case 'y':	# Dependent variable
	    call gt_sets (gt, gtlabel[axis], Memc[IC_LABELS(ic,2)])
	    call gt_sets (gt, gtunits[axis], Memc[IC_UNITS(ic,2)])
	    call amovr (y, z, npts)
	case 'f':	# Fitted values
	    call gt_sets (gt, gtlabel[axis], "fit")
	    call gt_sets (gt, gtunits[axis], Memc[IC_UNITS(ic,2)])
	    call icg_vector (ic, nl, x, z, npts)
	case 'r':	# Residuals
	    call gt_sets (gt, gtlabel[axis], "residuals")
	    call gt_sets (gt, gtunits[axis], Memc[IC_UNITS(ic,2)])
	    call icg_vector (ic, nl, x, z, npts)
	    call asubr (y, z, z, npts)
	case 'd':	# Ratio
	    call gt_sets (gt, gtlabel[axis], "ratio")
	    call gt_sets (gt, gtunits[axis], "")
	    call icg_vector (ic, nl, x, z, npts)
#	    iferr (call adiv$t (y, z, z, npts))
		call advzr (y, z, z, npts, icg_dvzr)
	case 'n':	# Linear component removed
	    call gt_sets (gt, gtlabel[axis], "non-linear component")
	    call gt_sets (gt, gtunits[axis], Memc[IC_UNITS(ic,2)])
	    xmin = IC_XMIN(ic)
	    xmax = IC_XMAX(ic)
	    a = nl_zeval (nl, real (xmin), 0.)
	    b = (nl_zeval (nl, real (xmax), 0.) - a) / (xmax - xmin)
	    do i = 1, npts
	        z[i] = y[i] - a - b * (x[i] - xmin)
	default:	# User axes types.
	    call malloc (label, SZ_LINE, TY_CHAR)
	    call malloc (units, SZ_LINE, TY_CHAR)
	    if (axis == 1) {
		call strcpy (Memc[IC_LABELS(ic,1)], Memc[label], SZ_LINE)
		call strcpy (Memc[IC_UNITS(ic,1)], Memc[units], SZ_LINE)
	        call amovr (x, z, npts)
	    } else {
		call strcpy (Memc[IC_LABELS(ic,2)], Memc[label], SZ_LINE)
		call strcpy (Memc[IC_UNITS(ic,2)], Memc[units], SZ_LINE)
	        call amovr (y, z, npts)
	    }
	    call icg_uaxes (axistype, nl, x, y, z, npts, Memc[label],
		Memc[units], SZ_LINE)
	    call gt_sets (gt, gtlabel[axis], Memc[label])
	    call gt_sets (gt, gtunits[axis], Memc[units])
	    call mfree (label, TY_CHAR)
	    call mfree (units, TY_CHAR)
	}
end


# ICG_DVZ -- Error action to take on zero division.

real procedure icg_dvzr (x)

real	x			# Numerator

begin
	return (1.)
end


# ICG_VECTOR -- Compute function values in linear or mag scale.

procedure icg_vector (ic, nl, x, z, npts)

pointer	ic
pointer	nl
real	x[npts]
real	z[npts]
int	npts

pointer	sp, dummy
int	i

real	log10()

begin
	call smark (sp)
	call salloc (dummy, npts, TY_REAL)
	call amovkr (0., Memr[dummy], npts)
	call nl_vector (nl, x, Memr[dummy], z, npts)
	call sfree (sp)

	if (IC_YAXIS(ic) == IC_MAG) {
	    do i = 1, npts
	        z[i] = IC_MAG0(ic) - 2.5 * log10 (z[i])
	}
end
