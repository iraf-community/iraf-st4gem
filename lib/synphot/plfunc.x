include "libsynphot.h"

# PLFUNC -- Produce a power law spectrum

procedure plfunc (refwave, index, units, nwave, wave, spec)

real	refwave		# i: reference wavelength (angstroms)
real	index		# i: power law exponent
char	units[ARB]	# i: physical units of power law
int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: output spectrum, in photlam
#--
int	iwave, done

string	negwave   "Reference wavelength cannot be negative"

string	badunits  "Unknown flux units for function pl"

int	is_magunit(), anytophot()
errchk	synphoterr

begin
	if (refwave <= 0.0)
	    call synphoterr (negwave, "pl")

	# Create power law spectrum

	if (is_magunit (units) == YES) {
	    do iwave = 1, nwave
		spec[iwave] = -2.5 * index * alog10 (wave[iwave] / refwave)
	} else {
	    do iwave = 1, nwave
		spec[iwave] = (wave[iwave] / refwave) ** index
	}

	# Convert to photlam units

	iferr {
	    done = anytophot (units, nwave, wave, spec)
	} then {
	    done = NO
	} else {
	    done = YES
	}

	if (done == NO)
	    call synphoterr (badunits, units)

end
