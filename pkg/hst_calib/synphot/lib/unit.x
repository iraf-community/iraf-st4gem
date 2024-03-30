include <synphot.h>

# UNIT -- compute the stimulus needed for unit response per cm^2 sec
#
# Computes the stimulus that is needed to produce a unit response 
# of one count per second per square centimeter.
# This is a measure of the inverse sensitivity of the passband.
#
# Input:
#	NWAVE		I4 number of wavelengths
#	WAVE(NWAVE)	R4 wavelengths
#	FILT(NWAVE)	R4 passband
#	FORM		C* determines the form of stimulus requested
# Output:
#	UNIT		R4 the unit stimulus.
#
# The unit stimulus is given in various representations by :
#
#   FLAM	U_lam = hc / [ Sum FILT(lam) lam * dlam ]
#   FNU		U_nu  = h  / [ Sum FILT(lam) dlam / lam ]
#   PHOTLAM	N_lam = 1  / [ Sum FILT(lam) dlam ]
#   PHOTNU	N_nu  = 1  / [ Sum FILT(lam) dlam * c / lam^2 ]
#   PHOTPIX	N_pix = 1  / [ Sum FILT(lam) ]
#   ABMAG	compute U_nu and convert to ABMAG
#   STMAG	compute U_lam and convert to STMAG
#   VEGAMAG	compute 1 / count rate on Vega
#   JY,MJY	compute U_nu and convert to JY,MJY
#
# Nov 1985 Keith Horne @ STScI - coded as STUNIT
# Feb 1987 KDH @ STScI - coded as UNIT with support for various forms.
# Mar 1988 KDH @ STScI - add VEGAMAG
# Dec 1988 E. Medeiros - SPP version
#
include <tbset.h>

real procedure unit ( nwave, wave, filt, form )

pointer	resample		# pointer to resampled vega flux in dynamic
				# memory

int    	nwave			# number of passband function samples
int	form_flag		# spectrum representation form flag
int  	nvega			# number of vega spectrum flux samples
int	i			# loop counter
int	formlist()

real	wave[ARB]		# passband wavelength set data
real	filt[ARB]		# passband transmisson data
real 	H 			# Planck's constant
real	C 			# speed of light in vacum   
real	vwave[MAXVEGA]		# vega flux wavelength set
real	vega[MAXVEGA]		# vega flux data
real	rate			# count rate function value
real	sumfilt			# integration function value

double	sum			# value of integration function

char	form[SZ_LINE]		# unit-response stimulus representation
char	vform[SZ_LINE]		# form of returned vega spectrum

bool	status			# execution status of specform (not used)

begin
 	
	# validate the input passband data
	if ( nwave <= 0 ){

	   call printf(" Number of wavelengths = 0\n")

	}else if ( wave[1] == wave[nwave] ) {

	   call printf(" Initial wavelength = final wavelength\n" )

	}

	#INITIALIZE local variables
	H = 6.62620E-27
	C = 2.997925E+18

	form_flag = formlist( form )

	   switch (form_flag) {
	      case 1,2: # spectrum representation is in PHOTLAM
                 unit = sumfilt ( nwave, wave, 0, filt )
                 if ( unit > 0.0 ) {
                    unit = 1. / unit
                 }
              case 3,4: # spectrum representation is in PHOTPIX
                 sum = 0.0D0
                 do i = 1, nwave {
	            if ( !IS_INDEFR (filt[i]) )
                       sum = sum + filt[i]
                 }
                 if ( sum > 0.0 ) {
                    unit = 1. / sum
                 }
              case 5,6: # spectrum representation is in FLAM
                 unit = sumfilt ( nwave, wave, 1, filt )
                 if ( unit > 0.0 ) {
	            unit = H * C / unit
                 }
              case 7,8: # spectrum representation is in FNU
                 unit = sumfilt ( nwave, wave, -1, filt )
                 if ( unit > 0.0 ) {
                    unit = H / unit
                 }
              case 9,10: # spectrum representation is in PHOTNU   
                 unit = sumfilt ( nwave, wave, -2, filt )
                 if ( unit > 0.0 ) {
                    unit = 1. / (C * unit)
                 }
              case 11,12: # spectrum representation is in JY     
                 unit = sumfilt ( nwave, wave, -1, filt )
                 if ( unit > 0.0 ) {
                    unit = H / unit
                    call specform ( 1, 1., unit, "FNU", unit, "JY", status )
                 }
              case 13,14: # spectrum representation is in MJY
                 unit = sumfilt ( nwave, wave, -1, filt )
                 if ( unit > 0.0 ) {
                    unit = H / unit
                    call specform ( 1, 1., unit, "FNU", unit, "MJY", status )
                 }
              case 15,16: # spectrum representation is in ABMAG
                 unit = sumfilt ( nwave, wave, -1, filt )
                 if ( unit > 0.0 ) {
                    unit = H / unit
                    call specform ( 1, 1., unit, "FNU", unit, "ABMAG", status )
                 }
              case 17,18: # spectrum representation is in STMAG
                 unit = sumfilt ( nwave, wave, 1, filt )
                 if ( unit > 0.0 ) {
                    unit = H * C / unit
                    call specform ( 1, 1., unit, "FLAM", unit, "STMAG", 
						status ) 
                 }
              case 19,20: # spectrum representation is in VEGAMAG
	          # Read in the vega spectrum.  rdvega returns actual 
	          # number of data points.
	          nvega = MAXVEGA
	          call rdvega( vwave, vega, nvega, vform )
                  call malloc ( resample, nwave, TY_REAL )
                  call linterp ( form, nvega, vwave, vega, nwave, wave, 
					                  Memr[resample] )
                  unit = rate ( nwave, wave, filt, Memr[resample], vform )
                  if ( unit > 0.0 ) {
                     unit = 2.5 * ALOG10( unit )
                  }
	          call mfree ( resample, TY_REAL ) 
	      case 21: # form is OBMAG
                 sum = 0.0D0
                    do i = 1, nwave {
	               if ( !IS_INDEFR (filt[i]) )
                          sum = sum + filt[i]
                    }
                    if ( sum > 0.0 ) {
                       unit = 1. / sum
                    }
                    call specform ( 1, 1., unit, "COUNTS", unit, "OBMAG", 
						status ) 

              default: # spectrum representation is in a unknown form
                  unit = 0.0
           }

        return
end
