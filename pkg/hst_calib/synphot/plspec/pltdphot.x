include	<gset.h>
include	"dphot.h"

#* HISTORY *
#* B.Simon	08-Jun-94	original
#* B.Simon	14-Jul-94	updated to handle indefs

# PLTDPHOT -- Plot photometric statistics and passbands

procedure pltdphot (gp, dphot)

pointer	gp		# i: Graphics descriptor
pointer	dphot		# i: Spectrum descriptor
#--
bool	plot
int	nwave, iband, istat
pointer	stim, pivot, fwhm, wave, band

begin
	if (dphot == NULL)
	    return

	# Plot photometric statistics

	stim = PHT_STIM(dphot)
	pivot = PHT_PIVOT(dphot)
	fwhm = PHT_FWHM(dphot)

	plot = false
	do istat = 1, PHT_NSTAT(dphot) {
	    if (! IS_INDEFR (Memr[stim])) {
		call gamove (gp, Memr[pivot]-0.5*Memr[fwhm], Memr[stim])
		call gadraw (gp, Memr[pivot]+0.5*Memr[fwhm], Memr[stim])
		call gmark (gp, Memr[pivot], Memr[stim], GM_CIRCLE, 2., 2.)

		plot = true
	    }

	    stim = stim + 1
	    pivot = pivot + 1
	    fwhm = fwhm + 1
	}

	# Plot the passbands

	if (plot) {
	    wave = PHT_WAVE(dphot)
	    band = PHT_BAND(dphot)
	    nwave = PHT_NWAVE(dphot)

	    do iband = 1, PHT_NBAND(dphot) {
		call gpline (gp, Memr[wave], Memr[band], nwave)
		band = band + nwave
	    }
	}

end
