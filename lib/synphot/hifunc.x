include "libsynphot.h"
define	HCK 		(H * C / KBOLTZ)

# HIFUNC -- Compute hydrogen absorption spectrum

procedure hifunc (temp, colden, nwave, wave, spec)

real	temp		# i: temperature, in degrees Kelvin
real	colden		# i: logarithm of column density (cm^-3)
int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: output spectrum, in photlam
#--
int	iwave, n, nn, n1, n2
real	z, zz, beta, x, q, qq, g1, g2, sum, gbf, xg, gff, a
real	ophyd, stimcor, tau, factor

string	negtemp  "Negative temperatures not allowed (Kelvin)"
string	negdens  "Negative column densities not allowed"

begin
	if (temp <= 0.0)
	    call synphoterr (negtemp, "hi")

	if (colden <= 0.0)
	    call synphoterr (negdens, "hi")

	# Compute black body spectrum 

	call bbfunc (temp, nwave, wave, spec)

	z = 1.0
	zz = z * z
	beta = 157890.0 * zz / temp

	# Correct black body for hydrogen absorption

	do iwave = 1, nwave {
	    if (spec[iwave] <= 0.0)
		next

	    # Compute hydrogen cross section

	    x = wave[iwave] * zz / 911.267

	    q = x ** (-0.3333333)
	    qq = q * q
	    g1 = 0.3458 * q
	    g2 = 0.03307 * qq

	    # Sum continua of first few levels explicitly

	    n1 = sqrt(x) + 1
	    n2 = n1 + 2
	    sum = 0.0
	    do n = n1, n2 {
		nn = n * n
		xg = x / nn - 0.5
		gbf = 1.0 + g1 * xg - g2 * (xg * xg + 1.25)
		sum = sum + exp (beta / nn) * gbf / (nn * n)
	    }

	    # Use continuum approximation for remaining levels

	    xg = x / beta
	    gff = 1. + g1 * (xg + 0.5) - 2.0 * g2 *(xg *(xg + 0.5) + 0.75)
	    a = n2 + 0.5
	    sum = sum + (exp(beta/(a*a)) - 1. + gff)/(beta + beta)

	    # Compute cross section and correct for stimulated emission

	    ophyd = 7.907e-18 * x * x * x * sum / zz * exp(-beta)
	    stimcor = 1.0 - exp (- HCK / (temp * wave[iwave]))
	    ophyd = ophyd * stimcor

	    if (ophyd <= 0.0) {
		spec[iwave] = 0.0

	    } else {
		# Compute optical depth from cross section

		if (colden < 80.0) {
		    tau = ophyd * 10.0 ** colden
		} else {
		    tau = ophyd * colden
		}

		# Compute radiative transfer

		if (tau < 1.0e-3) {
		    factor = tau * (1.0 - tau *(0.5 - tau / 6.0))
		} else {
		    factor = 1.0 - exp (- tau)
		}

		# Multiply black body spectrum by optical depth

		spec[iwave] = factor * spec[iwave]

	    }
	}

end
