include "../plspec/dspec.h"
include "../plspec/dsphot.h"
include "dratio.h"

#* HISTORY *
#* B.Simon	13-Jul-94	original

# STATRATIO -- Compute goodness of fit statistics for ratio of sphot to spec

procedure statratio (form, dspec, dsphot, dratio, chisq, bias, rms, nstat)

char	form[ARB]	# i: form of data
pointer	dspec		# i: spectrum descriptor
pointer	dsphot		# i: spectrophotmetry descriptor
pointer	dratio		# i: ratio descriptor
real	chisq		# o: chi squared
real	bias		# o: bias
real	rms		# o: root mean squared
int	nstat		# o: number of data points used in statistics
#--
int	nwave, snwave, snflux, pnwave, pnflux
pointer	sp, spec, sphot, perr, wave, swave, ospec, pwave, osphot, operr
real	const

int	is_magunit()

begin
	# Set stats to INDEF if we are not computing the ratio

	if (dratio == NULL) {
	    chisq = INDEFR
	    bias = INDEFR
	    rms = INDEFR
	    nstat = INDEFI

	    return
	}

	# Allocate memory for arrays interpolated on wavelength set

	wave = RAT_WAVE(dratio)
	nwave = RAT_NWAVE(dratio)

	call smark (sp)
	call salloc (spec, nwave, TY_REAL)
	call salloc (sphot, nwave, TY_REAL)
	call salloc (perr, nwave, TY_REAL)

	# Set constant used for extrapolated points according to type

	if (is_magunit (form) == NO) {
	    const = 0.0
	} else {
	    const = 100.0
	}

	# Read last spectrum

	snwave = SPC_NWAVE(dspec)
	snflux = SPC_NFLUX(dspec)
	swave = SPC_WAVE(dspec)
	ospec = SPC_FLUX(dspec) + snwave * (snflux - 1)

	# Interpolate spectrum on ratio's wavelength set

	call synextrap (const, snwave, Memr[swave], Memr[ospec], 
			nwave, Memr[wave], Memr[spec])

	# Read last spectrophotometry file and associated error

	pnwave = SPT_NWAVE(dsphot)
	pnflux = SPT_NFLUX(dsphot)
	pwave = SPT_WAVE(dsphot)
	osphot = SPT_FLUX(dsphot) + pnwave * (pnflux - 1)
	operr = SPT_ERROR(dsphot) + pnwave * (pnflux - 1)

	# Interpolate on ratio's wavelength set

	call synextrap (const, pnwave, Memr[pwave], Memr[osphot], 
			nwave, Memr[wave], Memr[sphot])

	if (IS_INDEFR(Memr[operr])) {
	    call amovkr (INDEFR, Memr[perr], nwave)
	} else {
	    call synextrap (const, pnwave, Memr[pwave], Memr[operr], 
			    nwave, Memr[wave], Memr[perr])
	}

	# Compute goodness of fit statistics

	call goodfit (Memr[sphot], Memr[perr], Memr[spec], nwave, form,
		      chisq, bias, rms, nstat)

	call sfree (sp)
end
