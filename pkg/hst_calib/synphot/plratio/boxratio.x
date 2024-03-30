include	<gset.h>
include	"../limit.h"

#* HISTORY *
#* B.Simon	13-Jul-94	Original

# BOXRATIO -- Plot and label the axes for the ratio plot

procedure boxratio (gp, inlimit, outlimit, form, title)

pointer	gp		# i: graphics descriptor
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# i: calculated plot limits
char	form[ARB]	# i: spectral form
char	title[ARB]	# i: plot title
#--
pointer	sp, ylabel
real	dif, temp, eqval, limit[4]

int	is_magunit()

string	xlabel    "WAVELENGTH (A)"
string	noextent  "No vertical extent in plot"

string	ratiolab  "Observed/Predicted"
string	diflabel  "Observed - Predicted (MAG)"

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (ylabel, SZ_FNAME, TY_CHAR)

	# Get ylabel from spectral form

	if (is_magunit (form) == NO) {
	    call strcpy (ratiolab, Memc[ylabel], SZ_FNAME)
	} else {
	    call strcpy (diflabel, Memc[ylabel], SZ_FNAME)
	}

	# Adjust plot limits

	call amovr (outlimit, limit, 4)

	dif = limit[TOP] - limit[BOTTOM]
	if (dif == 0.0)
	    call printerr_real (noextent, limit[BOTTOM])

	if (IS_INDEFR(inlimit[BOTTOM]))
	    limit[BOTTOM] = limit[BOTTOM] - 0.05 * dif

	if (IS_INDEFR(inlimit[TOP]))
	    limit[TOP] = limit[TOP] + 0.05 * dif

	# Swap limits for magnitude units

	if (is_magunit (form) == YES) {
	    temp = limit[TOP]
	    limit[TOP] = limit[BOTTOM]
	    limit[BOTTOM] = temp
	}

	# Set limits and plot axes

	call gswind (gp, limit[LEFT], limit[RIGHT], limit[BOTTOM], limit[TOP])
	call gseti (gp, G_TXQUALITY, GT_HIGH)

	call glabax (gp, title, xlabel, Memc[ylabel])

	# Plot zero line for spectra if form is not magnitude

	if (is_magunit (form) == NO) {
	    eqval = 1.0
	} else {
	    eqval = 0.0
	}

	call gamove (gp, limit[LEFT], eqval)
	call gadraw (gp, limit[RIGHT], eqval)

	# Unset cursor position for label

	call labelu
	call sfree (sp)
end
