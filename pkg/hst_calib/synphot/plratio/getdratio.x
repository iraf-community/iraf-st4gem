include	<tbset.h>
include "../plspec/dspec.h"
include "../plspec/dsphot.h"
include "dratio.h"

#* HISTORY *
#* B.Simon	12-Jul-94	original

# GETDRATIO -- Compute the ratio of spectrum and spectrophotometry files

procedure getdratio (wavetab, form, dspec, dsphot, dratio)

char	wavetab[ARB]	# i: wavelength table name
char	form[ARB]	# i: output form
pointer	dspec		# i: spectrum descriptor
pointer	dsphot		# i: spectrophotmetry descriptor
pointer	dratio		# o: ratio descriptor
#--
bool	logspace
int	mag, slen, plen, nwave, nflux, isphot, ispec
pointer	sp, tp, swave, pwave, wave, sphot, spec, flux, osphot, ospec
real	const, minwave, maxwave

data	logspace  / true /
real    setindef() # Set error cases to INDEF
extern	setindef
int	is_magunit(), tbtopn(), tbpsta()

begin
	# Check for existence of spectrum and spectrophotometry descriptors

	if (dspec == NULL || dsphot == NULL) {
	    dratio = NULL
	    return
	}

	# Check if this is a magnitude form

	mag = is_magunit (form)
	if (mag == NO) {
	    const = 0.0
	} else {
	    const = 100.0
	}

	# Get pointers to wavelength tables from descriptors

	swave = SPC_WAVE(dspec)
	slen = SPC_NWAVE(dspec)
	pwave = SPT_WAVE(dsphot)
	plen = SPT_NWAVE(dsphot)

	# Compute wavelength set for ratios

	if (wavetab[1] != EOS) {
	    tp = tbtopn (wavetab, READ_ONLY, NULL)
	    nwave = tbpsta (tp, TBL_NROWS)
	    call malloc (wave, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)

	} else {
	    minwave = max (Memr[swave], Memr[pwave])
	    maxwave = min (Memr[swave+slen-1], Memr[pwave+plen-1])

	    nwave = max (slen, plen)
	    call malloc (wave, nwave, TY_REAL)
	    call waveset (logspace, minwave, maxwave, nwave, Memr[wave])
	}

	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (sphot, nwave, TY_REAL)
	call salloc (spec, nwave, TY_REAL)

	# Allocate ratio descriptor

	nflux = SPC_NFLUX(dspec) * SPT_NFLUX(dsphot)

	call malloc (dratio, LEN_RATSTRUCT, TY_STRUCT)
	call malloc (flux, nflux*nwave, TY_REAL)

	RAT_NWAVE(dratio) = nwave
	RAT_NFLUX(dratio) = nflux
	RAT_WAVE(dratio) = wave
	RAT_FLUX(dratio) = flux

	# Compute ratios of spectrophotometic files to spectra

	osphot = SPT_FLUX(dsphot)
	do isphot = 1, SPT_NFLUX(dsphot) {
	    # Interpolate spectrophotometric file on wavelngth set

	    call synextrap (const, plen, Memr[pwave], Memr[osphot], 
			    nwave, Memr[wave], Memr[sphot])

	    ospec = SPC_FLUX(dspec)
	    do ispec = 1, SPC_NFLUX(dspec) {
		# Interpolate spectrum on wavelength set

		call synextrap (const, slen, Memr[swave], Memr[ospec], 
				nwave, Memr[wave], Memr[spec])

		# Compute difference if magnitude units, ratio if not

		if (mag == NO) {
		    call advzr (Memr[sphot], Memr[spec], Memr[flux], 
				nwave, setindef)
		} else {
		    call asubr (Memr[sphot], Memr[spec], Memr[flux], nwave)
		}

		ospec = ospec + slen
		flux = flux + nwave
	    }

	    osphot = osphot + plen
	}

	call sfree (sp)
end

# SETINDEF - Set error cases to INDEF

real procedure setindef (val)

real    val             # numerator in division
#--

begin
        return (INDEFR)
end
