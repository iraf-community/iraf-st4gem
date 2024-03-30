# SUMFILTSPEC - integrate passband[lambda] and spectrum[lambda] over
#		wavelength set
#
# SUMFILTSPEC = Sum [ FILT(I) SPEC(I) WAVE(I)**NPOW DWAVE(I) ]
#
#
# Input:
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	FILT(NWAVE)	R4 filter
#	SPEC(NWAVE)	R4 spectrum
#	NPOW		I4 wavelength exponent
#
# Feb 1987 Keith Horne @ STScI
# Dec 1988 E.Medeiros -- SPP version
# Oct 1990 Dave Bazell  add fix to return 0. if nwave <= 1
#
include <tbset.h>

real procedure sumfiltspec ( nwave, wave, npow, filt, spec )
int	nwave			# number of array elements
int	npow			# exponential power of wavelength
int	i			# loop counter
int	case_switch		# integration acceleration switch

real	wave[ARB]		# input wavelength set array
real	filt[ARB]		# input passband transmission array
real	spec[ARB]		# input spectrum flux array
real	w			# initial integral limit
real	wp			# terminal integral limit
real	wm			# trapaxoid midpoint

double	sum			# integral value of passband

pointer	sp, oldspec,oldfilt

begin

	# initialize the value of the integral
	sum = 0.0D0

	if ( nwave <= 1 ) 
	   return (0.)

	w = wave[1]
	wp = wave[2]

	if ( nwave > 1 ) { # a valid numeric function has been input

	   # indef values should make no contribution to the sum so set
	   # them to zero
	   call smark ( sp )
	   call salloc( oldspec, nwave, TY_REAL )
	   call salloc( oldfilt, nwave, TY_REAL )
	   do i = 1,nwave {
	      Memr[oldspec+i-1] = spec[i]
	      Memr[oldfilt+i-1] = filt[i]
	      if ( spec[i] >= INDEFR )
	         spec[i] = 0.
	      if ( filt[i] >= INDEFR )
	         filt[i] = 0.
	   }

	   # set case switch to reflect offset from npow value
	   case_switch = npow + 3

	   switch ( case_switch ) { # select case on power of wavelength

	      case 1: # npow is -2

	         # sum first function trapazoid
		 sum = filt[1] * spec[1] * (wp - w) / (w * w)

		 # sum next nwave -1 trapzoids
		 do i = 2, nwave-1 {
		    wm = w
		    w = wp
		    wp = wave[i+1]
		    sum = sum + filt[i] * spec[i] * (wp - wm) / (wp * wp)
		 }

		 # sum last trapzoid
		 sum = sum + filt[nwave] * spec[nwave] * (wp - w) / (wp * wp)

	      case 2: # npow is -1

		 # sum first trapazoid
		 sum = filt[1] * spec[1] * (wp - w) / w

		 # sum next nwave-1 trapazoids
		 do i = 2, nwave-1 {
		    wm = w
		    w = wp
		    wp = wave[i+1]
		    sum = sum + filt[i] * spec[i] * (wp - wm) / w
		 }

		 # sum last trapazoid
		 sum = sum + filt[nwave] * spec[nwave] * (wp - w) / wp

	      case 3: # npow is 0

		 # sum first trapazoid
		 sum = filt[1] * spec[1] * (wp - w)

		 # sum next nwave -1 trapazoids
		 do i = 2, nwave-1 {
		    wm = w
		    w =wp
		    wp = wave[i+1]
		    sum = sum + filt[i] * spec[i] * (wp - wm)
		 }

		 # sum last trapazoid
		 sum = sum + filt[nwave] * spec[nwave] * (wp - w)

              case 4: # npow is 1

		 # sum fisrt trapazoid
		 sum = sum + filt[1] * spec[1] * (wp - w) * w

		 # sum next nwave -1 trapazoids
		 do i = 2, nwave-1 {
		    wm = w
		    w = wp
		    wp = wave[i+1]
		    sum = sum +filt[i] * spec[i] * (wp - wm) * w
		 }

		 # sum last trapazoid
		 sum = sum + filt[nwave] * spec[nwave] * (wp - w) * wp

              case 5: # npow is 2

		 # sum first trapazoid
		 sum = filt[1] * spec[1] * (wp - w) * w * w
		
                 # sum next nwave-1 trapazoids
		 do i = 2, nwave-1 {
		    wm = w
		    w = wp
		    wp = wave[i+1]
		    sum = sum + filt[nwave] * spec[i] * (wp - w) * w * w
		 }

		 # sum last trapazoid
		 sum = sum + filt[nwave] * spec[nwave] * (wp - w) * wp * wp

              default:

		 # sum first trapazoid
		 sum = filt[1] * spec[1] * (wp - w) * (w ** npow)

		 # sum next nwave-1 trapazoids
		 do i = 2, nwave-1 {
		    wm = w
		    w = wp
		    wp = wave[i+1]
		    sum = sum + filt[i] * spec[i] * (wp - wm) * (wp ** npow)
		 }

		 # sum last trapazoid
		 sum = sum + filt[nwave] * spec[nwave] * (wp - w) * (wp ** npow)

	   }
	}

	# set integral value to relfect trapazoid sums
	sumfiltspec = 0.5D0 * sum

	do i = 1, nwave {
	   spec[i] = Memr[oldspec+i-1]
	   filt[i] = Memr[oldfilt+i-1]
	}

	call sfree( sp )

	return
end
