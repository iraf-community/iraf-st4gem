include	<mach.h>

# ONEINTERP -- Interpolate at a single wavelength in a spectrum

real procedure one_interp (idx, wavelen, nwave, wave, flux)

int	idx		# u: index of nearby point in spectrum
double	wavelen		# i: wavelength to interpolate at
int	nwave		# i: number of wavelengths in spectrum
real	wave[ARB]	# i: wavelength grid of spectrum
real	flux[ARB]	# i: flux of associated spectrum
#--
double	a, b, dwave
real	value

string	divzero  "Division by zero in oneinterp at point"
string	badvalue "Illegal data value in oneinterp at point"

begin
	# Find the index of the wavelength just above 
	# the interpolating wavelength

	if (idx > nwave)
	    idx = nwave

	while (idx > 1) {
	    if (wave[idx] < wavelen)
		break

	    idx = idx - 1
	}

	while (idx <= nwave) {
	    if (wave[idx] > wavelen)
		break

	    idx = idx + 1
	}

	# Return zero flux for wavelength outside of interpolation range

	if (idx == 1 || idx > nwave)
	    return (0.0)

	# Perform linear interpolation to get interpolated flux

	dwave = double(wave[idx]) - double(wave[idx-1])
	if (dwave == 0.0)
	    call syninterr (divzero, idx)

	a = (wave[idx] - wavelen) / dwave
	b = 1.0 - a

	value = a * flux[idx-1] + b * flux[idx]
	if (value > MAX_REAL)
	    call syninterr (badvalue, idx)

	return (value)
end
