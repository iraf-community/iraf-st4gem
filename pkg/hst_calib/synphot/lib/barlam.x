## BARLAM -- Computes bar wavelength of a broad bandpass
#
# barlam = exp( Int ln(lam) P(lam) dln(lam) / Int P(lam) dln(lam) )
#
# Input:
#	NWAVE	= Number of wavelengths
#	WAVE	= Wavelengths (positive,ascending)
#	P	= Thruputs (positive)
#
# Nov 1985 Keith Horne @ STScI
# Nov 1988 E. Medeiros -- SPP version
#
include <tbset.h>

real procedure barlam ( nwave, wave, p )

int	nwave				# number of throughput data points
real	wave[ARB]			# wavelength set
real	p[ARB]				# throughput data

int	i				# loop counter
int	im				# wavelength data pointer
int	ip				# throughput data pointer

double	wgt				# the data weight
double	sum				# partial sum
double  suml				# partial sum

begin

	if ( nwave <= 0. ) { # no data available set the fucntion value
			     #to reflect no evaluation
	   barlam = 0.0

	} else if ( wave[1] == wave[nwave] ) { # a nul o constant wavelength
					       # set function value to
	   barlam = wave[1]		       # reflect the wavelength value

	} else { # compute the center of the passband

	   # initialize the partial sums
	   sum = 0.D0
	   suml = 0.D0

	   do i = 1, nwave { # set pointers, wieghts, and running mean to
			     # reflect the weighted mean of the throughput
			     # data
	      im = max( i-1, 1)
	      ip = min( i+1, nwave )
	      wgt = p[i] * ( wave[ip] - wave[im] ) / ( wave[i] + wave[i] )
	      sum = sum + wgt
	      suml = suml + wgt * alog( wave[i] )
	   }
	   if ( sum > 0.D0 ) { # throughput data is not constant

	      barlam = Dexp( suml / sum )
	   
	   } else { # throughput data is constant

	      barlam = wave[1]
	   }	
	}
   	return
end
