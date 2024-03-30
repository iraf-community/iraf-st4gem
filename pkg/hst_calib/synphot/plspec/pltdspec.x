include	<gset.h>
include	"dspec.h"

# PLTDSPEC -- Plot a spectrum

procedure pltdspec (gp, form, dspec)

pointer	gp		# i: Graphics descriptor
char	form		# i: Output form
pointer	dspec		# i: Spectrum descriptor
#--
int	nwave, iflux
pointer	wave, flux

int	is_count()

begin
	if (dspec == NULL)
	    return

	# Plot  the spectra

	wave = SPC_WAVE(dspec)
	flux = SPC_FLUX(dspec)
	nwave = SPC_NWAVE(dspec)

	if (is_count (form) == NO) {
	    do iflux = 1, SPC_NFLUX(dspec) {
		call gpline (gp, Memr[wave], Memr[flux], nwave)
		flux = flux + nwave
	    }

	} else {
	    do iflux = 1, SPC_NFLUX(dspec) {
		call histplot (gp, Memr[wave], Memr[flux], nwave)
		flux = flux + nwave
	    }
	}

end
