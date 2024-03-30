## FWHMLAM -- gaussian equivalent FWHM of passband
#
# Computes FWHM of gaussian passband with same rms width as FILT.
#
# Input:
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	FILT(NWAVE)	R4 passband
# Output:
#	FWHMLAM		R4 equivalent FWHM bandwidth
#
# Feb 1987 KDH @ STScI
# Nov 1988 E.Medeiros -- SPP version
#
include <tbset.h>

real procedure fwhmlam ( nwave, wave, filt )

int	nwave			# number of array elements

real	wave[ARB]		# input wavelength set array
real	filt[ARB]		# input passband transmission array

real	rmslam			# root mean square average wavelength

begin

	#set the full width half maximum to reflect the product of the
	#root mean square average wavelength and the 1 sigma range
	fwhmlam = SQRT( 8. * alog(2.)) * rmslam( nwave, wave, filt )
	return

end
