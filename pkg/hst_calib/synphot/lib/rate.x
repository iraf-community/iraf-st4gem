## RATE -- Computes the count rate produced by given spectrum in given passband.
# 	  output is in counts / second / cm^2
#
# Input:
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	FILT(NWAVE)	R4 passband
#	SPEC(NWAVE)	R4 spectrum
#	FORM		C* form of input spectrum data
#
# Output:
#	RATE 		R4 count rate (counts/s/cm^2)
#
# The count rate is given for various spectrum representations by :
#
#   FLAM	RATE	= [ Sum FILT(lam) * SPEC(lam) * dlam * lam ] / (hc)
#   FNU			= [ Sum FILT(lam) * SPEC(lam) * dlam / lam ] / h
#   PHOTLAM		= [ Sum FILT(lam) * SPEC(lam) * dlam ]
#   PHOTNU		= [ Sum FILT(lam) * SPEC(lam) * dlam / lam^2 ] * c
#   PHOTPIX		= [ Sum FILT(lam) * SPEC(lam) ]
#   ABMAG	convert to FNU then calculate as above
#   STMAG	convert to FLAM then calculate as above
#   VEGAMAG	convert to FLAM then calculate as above
#   OBMAG	convert to PHOTPIX then calculate as above
#
# Nov 1985 Keith Horne @ STScI - coded as STUNIT
# Feb 1987 Keiht Horne @ STScI - coded as UNIT with support for various forms.
# Dec 1988 E.Medeiros - SPP version
# Nov 1990 Dave Bazell  Change VEGAMAG conversion from FNU to FLAM and back
#
include	<tbset.h>

real procedure rate ( nwave, wave, filt, spec, form )

int	nwave		   	# number of array elements
int	form_flag		# spectrum representation form flag
int	i			# loop counter
int	formlist()
                            
real	wave[ARB] 		# wavelength set array
real	filt[ARB]		# passband(throughput) data array
real	spec[ARB]		# input flux data array
real 	H 			# Planck's constant
real	C 			# speed of light in vacum
real	sumfiltspec		# integration function value
real	sum			# local integration value

char	form[SZ_LINE]		# spectrum representation form

bool	status			# execution status of procedure rate
				#	success		=> 0
				#	failure		=> != 0

begin

	#INITIALIZE local variables
	H = 6.62620E-27
	C = 2.997925E+18
	status	=	true

	form_flag = formlist ( form )

	   switch (form_flag) {
	      case 1,2: # spectrum representation is in PHOTLAM

                 rate = sumfiltspec ( nwave, wave, 0, filt, spec )

              case 3,4: # spectrum representation is in PHOTPIX or COUNTS 

                 sum = 0.D0
                 do i = 1, nwave {
	            if ( !IS_INDEFR (filt[i]) && !IS_INDEFR (spec[i]) )
                       sum = sum + filt[i] * spec[i]
                 }
                 rate = sum

              case 5,6: # spectrum representation is in FLAM

                 rate = sumfiltspec ( nwave, wave, 1, filt, spec )
                 rate = rate / ( H * C )

              case 7,8: # spectrum representation is in FNU

                 rate = sumfiltspec ( nwave, wave, -1, filt, spec )
                 rate = rate / H

              case 9,10: # spectrum representation is in PHOTNU   

                 rate = sumfiltspec ( nwave, wave, -2, filt, spec )
                 rate = rate * C

              case 11,12: # spectrum representation is in JY

                 rate = sumfiltspec ( nwave, wave, -1, filt, spec )
                 rate = rate / H
                 call specform ( 1, 1., rate, "jy", rate, "fnu", status )

              case 13,14: # spectrum representation is in MJY

                 rate = sumfiltspec ( nwave, wave, -1, filt, spec )
                 rate = rate / H
                 call specform ( 1, 1., rate, "mjy", rate, "fnu", status )

              case 15,16: # spectrum representation is in ABMAG

                 call specform ( nwave, wave, spec, "abmag", spec, "fnu",
							status )
                 rate = sumfiltspec ( nwave, wave, -1, filt, spec )
                 rate = rate / H
                 call specform ( nwave, wave, spec, "fnu", spec, "abmag",
							status )

              case 17,18: # spectrum representation is in STMAG

                 call specform ( nwave, wave, spec, "stmag", spec, "flam",
							status )
                 rate = sumfiltspec ( nwave, wave, 1, filt, spec )
                 rate = rate / (H * C)
                 call specform ( nwave, wave, spec, "flam", spec, "stmag",
							status )

              case 19,20: # spectrum representation is in VEGAMAG

                 call specform ( nwave, wave, spec, "vegamag",
					spec, "flam", status )
                 rate = sumfiltspec ( nwave, wave, 1, filt, spec )
                 rate = rate / (H * C)
                 call specform ( nwave, wave, spec, "flam",
				 spec, "vegamag", status )

	      case 21,22: # spectrum representation is in OBMAG
	         call specform( nwave, wave, spec, "obmag", spec,
	                        "counts", status )
                 sum = 0.D0
                 do i = 1, nwave {
	            if ( !IS_INDEFR (filt[i]) && !IS_INDEFR (spec[i]) )
                       sum = sum + filt[i] * spec[i]
                 }
	         call specform( nwave, wave, spec, "counts", spec,
	                        "obmag", status )
                 rate = sum


              default: # spectrum representation is in a unknown form
	         
                 status = false
	         rate = 0.0
           }
        return
end
