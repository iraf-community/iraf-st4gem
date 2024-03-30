include "libsynphot.h"

# BBFUNC -- Compute a black body spectrum

procedure bbfunc (temp, nwave, wave, spec)

real	temp		# i: black body temperature, in Kelvin
int	nwave		# i: length of wavelength and spectrum arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: black body spectrum, in photlam units
#--
int	iwave
real	factor, x, c1, c2

data 	c1 / 1.43883e8 /	# hc/k [Angstroms-Kelvin]
data 	c2 / 1.95722e5 /	# numerical factor to make units come out OK

begin
	do iwave = 1, nwave {
	    x = wave[iwave] * temp

	    if (x <= 0.0) {
		spec[iwave] = 0.0

	    } else {
		# Black body spectrum in cgs units

		x = c1 / x
		if (x < 1.0e-4) {
		    factor = 2.0 / (x * (x + 2.0))
		} else if (x < 85.0) {
		    factor = 1.0 / (exp (x) - 1.0)
		} else {
		    factor = 0.0
		}

		x = x * temp / c2
		x = factor * x * x * x

		# Convert from cgs to photlam

		spec[iwave] = x / (H * wave[iwave])
	    }
	}
end

