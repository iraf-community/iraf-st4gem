# EFFSTIM -- Computes effective stimulus of passband and spectrum
#
# Computes the effective stimulus of a given spectrum in given passband.
# The effective stimulus is the count rate times the unit stimulus.
#
# Input:
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	FILT(NWAVE)	R4 passband
#	SPEC(NWAVE)	R4 spectrum
#	FORM		C* form of input spectrum data
#
# Output:
#	EFFSTIM		R4 effective stimulus in the representation form
#
# The effective stimulus is given for various spectrum representations by :
#
#   FLAM	
#   FNU			
#   PHOTLAM		
#   PHOTNU		
#   PHOTPIX		
#   ABMAG	
#   STMAG	
#   VEGAMAG
#
# Feb 1987 Keith Horne @ STScI 
# Dec 1988 E.Medeiros - SPP version
# Nov 1990 Dave Bazell  Change forms to lower case rather than upper
# Nov 1993 B.Simon	Call inisyntab,clssyntab

include	<tbset.h>

real procedure effstim ( nwave, wave, filt, spec, form )

int	nwave		    	# number of array elements
int	strmatch		# string match function value
                            
real	wave[ARB] 		# wavelength set array
real	filt[ARB]		# passband(throughput) data array
real	spec[ARB]		# input flux data array
real	rate			# expected count rate function value
real	unit			# effective unit stimulus function value
real 	count_rate		# expected count rate 
real	result			# effective stimulus

char	form[SZ_LINE]		# spectrum representation form
char	mag_pat[SZ_LINE,2]	# magnitude representation pattern

begin

	call inisyntab

	#initialize the magnitude pattern string
	call strcpy ( "MAG", mag_pat[1,1], 3 )
	call strcpy ( "mag", mag_pat[1,2], 3 )

	if ( nwave <= 0 ) { # wavelength data is present

	   call printf(" Number of wavelengths = 0\n")
	   result = 0.0

	}else if ( wave[1] == wave[nwave] ) { # wavelength data is not constant

	   call printf(" Initial wavelength = final wavelength\n" )
	   result = 0.0

	}else{ # input data has been validated

	   # calculate the expected count rate
	   count_rate = rate ( nwave, wave, filt, spec, form )

	   if ( count_rate <= 0.0 ) { # spectrum data is negative
	      result = 0.0

	   } else {

	      result = unit ( nwave, wave, filt, form )

	      if ( result <= 0.0 ) { # passband transmission is negative
		 result = 0.0

	      } else {
		 # convert the spectrum data representation to lower case
	         call strlwr ( form )

	         if ( strmatch ( form, mag_pat[1,1] ) > 0 ||
                      strmatch ( form, mag_pat[1,2] ) > 0 ) { # a magnitude
			# spectrum representation has been detected

	            result = result - 2.5 * alog10( count_rate )

                 }else{ # a spectrum flux representation detected

                    result = result * count_rate

                 }
              }
           }
	}

	call clssyntab
	return (result)
end
