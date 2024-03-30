include	<gset.h>
include	"dband.h"

# PLTDBAND -- Plot a set of passbands

procedure pltdband (gp, ylog, dband)

pointer	gp		# i: Graphics descriptor
bool	ylog		# i: Plot on logarithmic scale?
pointer	dband		# i: Spectrum descriptor
#--
int	nwave, iband
pointer	sp, wave, band, band2

begin
	if (dband == NULL)
	    return

	# Get pointers to arrays from structure

	wave = BND_WAVE(dband)
	band = BND_BAND(dband)
	nwave = BND_NWAVE(dband)

	call smark (sp)
	call salloc (band2, nwave, TY_REAL)

	do iband = 1, BND_NBAND(dband) {
	    # If this is a log plot, all non-positive values 
	    # must be set to INDEF

	    if (ylog) {
		call choplim (0.0, INDEFR, Memr[band], Memr[band2], nwave)
		call gpline (gp, Memr[wave], Memr[band2], nwave)

	    } else {
		call gpline (gp, Memr[wave], Memr[band], nwave)
	    }

	    band = band + nwave
	}

	call sfree (sp)
end
