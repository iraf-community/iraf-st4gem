# RATESPEC -- compute pixellated count-rate spectrum
#
# Computes count rate for given spectrum in `spectrograph' pixels.
#
# The passband of a spectrograph pixel is specified in one of two ways:
# If WAVE1(I) < WAVE2(I), the passband is given by PASS with sharp edges
#		at WAVE1 and WAVE2.
# If WAVE1(I) > WAVE2(I), the passband is gaussian with center WAVE1 and
#		full-width-half-maximum WAVE2.
#
# Output spectrum is counts/second/cm^2 .
#
# Input:
#	NPIX		I4 number of pixels
#	WAVE1(NPIX)	R4 pivot wavelength or lower wavelength limit
#	WAVE2(NPIX)	R4 fwhm wavelength or upper wavelength limit
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	PASS(NWAVE)	R4 passband
#	SPEC(NWAVE)	R4 spectrum
#	FORM		C* form of input spectrum data
#
# Output:
#	PIXRATE(NPIX)	R4 count rate (counts/s/cm^2)
#--
# Feb 1987 Keith Horne @ STScI - first version implements sharp pixel
#		boundaries but not gaussian passbands.
# Dec 1988 E. Medeiros	SPP version
include	<tbset.h>

procedure ratespec ( npix, wave1, wave2, pixrate, nwave, wave, pass,
				spec, form )

int	npix		# number of pixels
int	nwave		# number of passband and spectrum flux values
int	ipix		# pixel geometry loop counter
int	i		# loop counter

real	wave1[ARB]	# pivot wavelength or lower wavelength limit
real	wave2[ARB]	# fwhm wavelength or upper wavelength limit
real	pixrate[ARB]	# count rate per pixel ( counts/s/cm^2 )
real	wave[ARB]	# wavelength samples on passband and spectrum
real	pass[ARB]	# specified passband transmission
real	spec[ARB]	# specified spectrum fluxes            
real	pixwave[2]	# pixel boundary wavelength array
real	pixpass[2]	# pixel boundary passband transmission array
real	pixspec[2]	# pixel boundary flux array
real	rate		# count rate function value

char	form[SZ_LINE]	# the spectrum flux data representation

bool	warning		# flag to indicate sending of non-gaussion
			# warning message

begin

	# validate input data 
	if ( nwave <= 0 ) {

	   call printf("  No input passband data\n")

	}else if ( wave[1] == wave[nwave] ) {

	   call printf(" No input passband data\n")

	}else if ( npix <= 0 ) {

	   call printf(" No detector pixels\n")

	}else { # proceed through pixel geometry computing the count
		# rate for each pixel of detector

	   # set non-gaussian pixel geometry request
	   warning = true

	   do ipix = 1, npix { # loop through pixel geometry arrays

	      # get wavelengths at pixel boundaries
	      if ( wave1[ipix] > wave2[ipix] ) { # gaussian pixel
					#geometry requested
	         # spoof edges of pixel boundaries
	         pixwave[1] = wave1[ipix] - 0.5 * ABS( wave2[ipix] )
	         pixwave[2] = wave1[ipix] + 0.5 * ABS( wave2[ipix] )

	         if ( warning ) { # non-gaussian pixel message

	            call printf(" WARNING: gaussian pixel passbands not yet\n")
		    call printf(" implemented. Using appropriate boxcar\n")
		    call printf(" boundaries instead\n")

		    warning	=	false
		 }
              }else{ # set the boundary wavelengths to boxcar geometry
		 pixwave[1] = wave1[ipix]
		 pixwave[2] = wave2[ipix]
	      }
	      # locate pixel boundaries in the wavelength and passband
	      # arrays
	      call linterp ("", nwave, wave, pass, 1, pixwave[1], pixpass[1])
	      call linterp ("", nwave, wave, pass, 1, pixwave[2], pixpass[2])
	      call linterp (form, nwave, wave, spec, 1, pixwave[1], pixspec[1])
	      call linterp (form, nwave, wave, spec, 1, pixwave[2], pixspec[2])

	      # compute pixel count rate using function rate
      	      pixrate[ipix] = rate ( 2, pixwave, pixpass, pixspec, form )
           }
	   return
        } # invalide input, return nul pixrate set
        do i =1, npix {
	   pixrate[i] = 0.0
        }
	return
end
