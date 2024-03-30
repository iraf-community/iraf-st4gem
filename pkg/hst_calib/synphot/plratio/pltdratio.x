include	<gset.h>
include	"dratio.h"

#* HISTORY *
#* B.Simon	12-Jul-94	original

# PLTDRATIO -- Plot a spectrum

procedure pltdratio (gp, form, dratio)

pointer	gp		# i: Graphics descriptor
char	form		# i: Output form
pointer	dratio		# i: Spectrum descriptor
#--
int	nwave, iflux
pointer	wave, flux

int	is_count()

begin
	if (dratio == NULL)
	    return

	wave = RAT_WAVE(dratio)
	flux = RAT_FLUX(dratio)
	nwave = RAT_NWAVE(dratio)

	if (is_count (form) == NO) {
	    do iflux = 1, RAT_NFLUX(dratio) {
		call gpline (gp, Memr[wave], Memr[flux], nwave)
		flux = flux + nwave
	    }

	} else {
	    do iflux = 1, RAT_NFLUX(dratio) {
		call histplot (gp, Memr[wave], Memr[flux], nwave)
		flux = flux + nwave
	    }
	}

end
