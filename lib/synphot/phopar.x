include	<synphot.h>

# PHOPAR -- Calculate the four photometric parameters for HST

procedure phopar (nwave, wave, thruput, phot)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: throughput array
real	phot[4]		# o: photometric parameters
#--
real	funit(), pivlam(), rmslam()

begin
	# The four photometric parameters are:
	# [1] the inverse sensitivity
	# [2] the zero point parameter
	# [3] the pivot wavelength
	# [4] the rms bandwidth

	phot[1] = funit (HSTAREA, nwave, wave, thruput)
	phot[2] = STZERO
	phot[3] = pivlam (nwave, wave, thruput)
	phot[4] = rmslam (nwave, wave, thruput)
end
