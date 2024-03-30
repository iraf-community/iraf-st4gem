## AVGLAM -- Computes average wavelength of a passband
#
# Computes average wavelength of a broad bandpass
# avglam = Int lam P(lam) dlam / Int P(lam) dlam
#
# Input:
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	FILT(NWAVE)	R4 passband
# Output:
#	AVGLAM		R4 average wavelength
#
# Nov 1985 Keith Horne @ STScI
# Feb 1987 KDH @ STSCI - new version using SUMFILT
# Nov 1988 E.Medeiros -- SPP version
#
include <tbset.h>

real procedure avglam ( nwave, wave, filt )
 
int	nwave			# number of array elements

real	wave[ARB]		# input wavelength set array
real	filt[ARB]		# input passband transmission array
real	sum			# integral value of passband
real 	sum2			# integral value squared of passband
real	sumfilt			# numerical integration function value

begin

	if ( nwave <= 0 ) { # no data available set the function value to 
			    # reflect no evaluation
	   avglam = 0.0

	} else if ( wave[1] == wave[nwave] ) { # a nul or constant wavelength 
					       # set set function value to 
	   avglam = wave[1]		       # reflect the wavelength value

	} else { # compute numerical passband integral using SUMFILT with
		 # NPOW = 0 ( the general case )
           
	   sum = sumfilt ( nwave, wave, 0, filt )

 	   if ( sum <= 0.0 ) { # nul passband integral detected set the average
			       # wavelength to reflect the any wavelength value
	      avglam = wave[1]

	   } else { # compute numerical passband integral using SUMFILT with
		    # NPOW = 1

	      sum2 = sumfilt ( nwave, wave, 1, filt )
	      avglam =  sum2 / sum 

	   }
	}
	return
end
