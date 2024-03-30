include "libsynphot.h"

# BBSI -- Compute a black body spectrum

procedure bbsi (temp, nwave, wave, spec)

real    temp        # i: black body temperature, in kelvin
int     nwave       # i: length of wavelength and spectrum arrays
real    wave[ARB]   # i: wavelengths at which spectrum is computed (Angstroms)
real    spec[ARB]   # o: black body spectrum in photlam per square arcsec
#--
int	iwave
real	x
real	lam
real 	c
real 	h
real 	k
real 	f
real	llam_SI()

#	>>> af = 0.01 * 0.01	# per m^2  -->  per cm^2
#	>>> af
#	0.0001
#	>>> sf = 206265.0 * 206265.0
#	>>> sf = 1/sf
#	>>> sf			# per sr  -->  per sqarcsec
#	2.3504386381829067e-11
#	>>> af * sf
#	2.3504386381829069e-15
#	>>> af * sf * 1.0e-10	# per m  -->  per Angstrom
#	2.3504386381829069e-25
#
# per m^2 per steradian per m -->
# per cm^2 per square arcsecond per A
data	f / 2.3504386381829069e-25 /	


begin
	#* si.c 		anands@stsci.edu*/
	c	= 2.997925e+8
	h	= 6.6262e-34
	k	= 1.38062e-23

	do iwave = 1, nwave {
	    x = wave[iwave]

	    if (x <= 0.0) {
		spec[iwave] = 0.0

	    } else {
		# Black body spectrum in 
		# number of photons  per square cm  per second
		# per square arc second  per angstrom

		lam = x * 1.0E-10 # convert to meters

		spec[iwave] = f * llam_SI(lam, temp) / (h * c / lam)
	    }
	}
end



# derivative of nicmos cgi etc C code blackbody.c      anands@stsci.edu

real procedure llam_SI(lam, T)
	#* From Kaye and Laby 14th Ed 1973 p. 76 */
	#* returns Power area-1 solidangle-1 wavelength_range-1 */
	#* in strictly SI units

real lam          # i: wavelength in meters
real T            # i: temp in kelvin

real exponent
real c1prime
real c2
real c
real h
real k
real result

define TOO_SMALL  1.0E-4
define TOO_BIG    85.0

begin
	#* si.c 		anands@stsci.edu*/
	c	= 2.997925e+8
	h	= 6.6262e-34
	k	= 1.38062e-23

	c1prime = 2 * h * c * c	#* Power * unit area / steradian */
	c2 = h * c / k

# Some sleight of hand here is required in order to avoid
# floating point overflows in extreme cases. This logic, and
# the associated magic numbers, were taken from bbfunc.
	exponent = c2 / (lam * T)

	if (exponent < TOO_SMALL) {
            result = (2 * c1prime * (lam**-5.0)) / (exponent * (exponent + 2.0) )
	} else if (exponent < TOO_BIG) {
            result = (c1prime * (lam**-5.0)) / (exp(exponent) - 1.0)
	} else {
	    result = 0.0
	}

#Original code without checks.
#	expfactor = exp(exponent)
#	result =  (c1prime * (lam**-5.0)) / (expfactor - 1.0)
#...........................................................
#	call printf("\n Lam, result: %g %g ")
#	call pargr(lam)
#	call pargr(result)
#	call flush(STDOUT)

	return (result)
end
