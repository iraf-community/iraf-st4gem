## EXTMAG -- Evaluation of Seaton's extinction function 
#
# seaton's paper in m.n.r.a.s. vol 187, page 75p (1979).  
# the formulae are based on an adopted value of R = 3.20.
#
# note that seaton's representation of of the interstellar reddening law
# differs substantially from schild's representation (astron. j. 82, 339,
# table ii, 1977) in the region of overlap.  schild also adopted r = 3.20.
#
# for wavelengths > 3704 angstroms, the program interpolates 
# linearly in 1/lambda c in seaton's table 3.  
# for wavelengths < 3704 angstroms, the program uses the formulae 
# from seaton's table 2.  
# the formulae match at the endpoints of their respective intervals. 
# there is a mismatch of 0.009 mag/ebmv at nu=2.7 (lambda=3704 angstroms).
# seaton's tabulated value of 1.44 mags at 1/lambda = 1.1 may be in error;
# 1.64 seems more consistent with his other values.
#      
# wavelength range allowed is 0.1 to 1.0 microns.
# calls to the subroutine for wavelengths outside this range
# result in extrapolated values for the extinction being returned.
#
# sources:
#	lambda < 1000		same as lambda = 1000.
#	1000 < lambda < 3704	Seaton(1979) MNRAS 187,73p.
#	3704 < lambda < 10,000	Nandy(1975) A+A 44, 195. (corrected to R=3.2)
#	10000 < lambda		extrapolate linearly in 1/lam (can be improved)
#
# Input:
#	WAVE	= wavelength (Angstroms)
#	EBMV	= extinction parameter (mags) (can be negative)
#
# Output:
#	EXTMAG	= extinction in magnitudes
#
# Mar 1986 Keith Horne @ STScI -- adapted from R.WADE routine SEATON
# Dec 1988 E. Medeiros -- SPP version
# Dec 1989 Dave Bazell; correct error in etable, remove restriction at
# 	short wavelengths, add linear extrapolation at long wavelengths

include <tbset.h>

real procedure extmag ( wave, ebmv )

int i				# loop pointer
int ntable			# number of entrys in Seaton's table
data ntable/19/

real wave			# wavelength for Seaton's evaluation
real ebmv			# extinction parameter in Vmag
real x				# wavelength in microns truncated on 1000A

real xtable[19]			# tabulated inverse wavelengths
real etable[19]			# tabulated extinctin at E(B-V)=1.

data (xtable(i),i=1,19)/0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
			1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
			2.2, 2.3, 2.4, 2.5, 2.6, 2.7/

data (etable(i),i=1,19)/0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44,
     			2.66, 2.88, 3.14, 3.36, 3.56, 3.77,
     			3.96, 4.15, 4.26, 4.40, 4.52, 4.64/

begin

	if ( wave <= 0. )
	   return

	# convert wave in angstroms to 1/microns
	x = 10000.0 / wave

	# Infrared - extend optical results linearly to 0 at 1/lam = 0
	if ( x <= 1.0)
	   extmag = etable[2] * x * x

	else if ( x < 2.7 ) { # optical region interpolates
		 	  # in Seaton's table 3

	   # extract extmag by linear interpolation through Seaton's table
	   call syninterp ( ntable, xtable, etable, 1, x, extmag )

	}else if ( x < 3.65 ) {# ultraviolet uses analytic formulae from
			       # Seaton's table 2

	   extmag = 1.56 + 1.048 * x + 1.01 /
			(( x - 4.6 ) * ( x - 4.6 ) + 0.280)

	}else if ( x < 7.14 ) { # more ultraviolet

	   extmag = 2.29 + 0.848 * x + 1.01 /
			(( x - 4.6 ) * ( x - 4.6 ) + 0.280)

	}else if ( x <= 10. ) { # and more untraviolet still

	   extmag = 16.17 + x * ( -3.20 + 0.2975 * x )

	} else {

	   x = min( x, 50.)
	   extmag = 16.17 + x * ( -3.20 + 0.2975 * x )
	}

	# rescale the magnitude extinction to the Vmag extinction
	extmag = extmag * ebmv

	return
end
