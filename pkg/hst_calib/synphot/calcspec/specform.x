include	<synphot.h>
define	TWONFORM 22

## SPECFORM -- Convert from one data form to another.
#
# Input:
#	NWAVE	= I4 number of wavelengths
#	WAVE	= R4 pivot wavelengths (in Angstroms)
#	SPECIN	= R4 input spectrum data
#	INFORM	= C* form of spectrum data on input (e.g. 'FNU')
#	OUTFORM	= C* form of spectrum desired on output
# Output:
#	SPECOUT	= R4 output spectrum data
#	OUTFORM	= C* form of output spectrum data 
#
# Spectra may be converted among any of the following forms:
# Flux densities:
#	FLAM	= ergs cm-2 s-1 A-1
#	FNU	= ergs cm-2 s-1 Hz-1
#	JY	= 10^-23 ergs cm-2 s-1 Hz-1
#	MJY	= 10^-26 ergs cm-2 s-1 Hz-1
#	PHOTNU	= photons cm-2 s-1 Hz-1
#	PHOTLAM = photons cm-2 s-1 A-1
#	COUNTS  = photons cm-2 s-1 pixel-1 (requires monotonic wavelengths)
# Magnitude systems:
#	STMAG	= -2.5 log( FLAM ) - 21.10
#	ABMAG	= -2.5 log( FNU ) - 48.60
#	VEGAMAG	= -2.5 log( f / f(vega) )
#	OBMAG	= -2.5 log( COUNTS )
# Jan 1987 Keith Horne @ STScI - original version
# E.Medeiros - SPP version
# Dec. 1989 Dave Bazell - Change PHOTPIX to COUNTS and add OBMAG
#
include	<tbset.h>

procedure specform ( nwave, wave, specin, inform, specout, outform, status )

pointer	vdata			#pointer to start of dynamic memory
				#containing converted Vega data

int	nwave			#number of array elements
int	nvega			#number of vega data array elements
int	inform_flag		#indicator switch for input data form
int	outform_flag		#indicator switch for output data form
int	im			#array element pointer
int	ip			#array element pointer
int	i			#loop pointer

real	wave[ARB] 		#wavelength set array
real	specin[ARB]		#input flux data array
real	specout[ARB]		#output flux data array
real	factor			#convertion factor
real 	H 			#Planck's constant
real	C 			#speed of light in vacum
real	vwave[291]		#vega flux wavelength set
real	vega[291]		#vega flux data

char	inform[SZ_LINE]		#input flux form
char	outform[SZ_LINE]	#output flux form
char	forms[SZ_LINE,TWONFORM] #flux forms

bool	status			#execution status word
				# success	=>	true
				# failure	=>	false
bool	streq			#boolian string match function value 

 data nvega/291/

 data (vwave(i),i=1,32)/1100.,3300.,3325.,3350.,3375.,3400.,3425.,3450.,
             3475.,3500.,3525.,3550.,3575.,3600.,3625.,3650.,
             3675.,3700.,3725.,3750.,3775.,3800.,3825.,3850.,
             3875.,3900.,3925.,3950.,3975.,4000.,4025.,4050./

 data (vwave(i),i=33,64)/4075.,4100.,4125.,4150.,4175.,4200.,4225.,4250.,
             4275.,4300.,4325.,4350.,4375.,4400.,4425.,4450.,
             4475.,4500.,4525.,4550.,4575.,4600.,4625.,4650.,
             4675.,4700.,4725.,4750.,4775.,4800.,4825.,4850./

 data (vwave(i),i=65,96)/4875.,4900.,4925.,4950.,4975.,5000.,5025.,5050.,
             5075.,5100.,5125.,5150.,5175.,5200.,5225.,5250.,
             5275.,5300.,5325.,5350.,5375.,5400.,5425.,5450.,
             5475.,5500.,5525.,5550.,5575.,5600.,5625.,5650./

 data (vwave(i),i=97,128)/5675.,5700.,5725.,5750.,5775.,5800.,5825.,5850.,
             5875.,5900.,5925.,5950.,5975.,6000.,6025.,6050.,
             6075.,6100.,6125.,6150.,6175.,6200.,6225.,6250.,
             6275.,6300.,6325.,6350.,6375.,6400.,6425.,6450./

 data (vwave(i),i=129,160)/6475.,6500.,6525.,6550.,6575.,6600.,6625.,6650.,
             6675.,6700.,6725.,6750.,6775.,6800.,6825.,6850.,
             6875.,6900.,6925.,6950.,6975.,7000.,7025.,7050.,
             7075.,7100.,7125.,7150.,7175.,7200.,7225.,7250./

 data (vwave(i),i=161,192)/7275.,7300.,7325.,7350.,7375.,7400.,7425.,7450.,
             7475.,7500.,7525.,7550.,7575.,7600.,7625.,7650.,
             7675.,7700.,7725.,7750.,7775.,7800.,7825.,7850.,
             7875.,7900.,7925.,7950.,7975.,8000.,8025.,8050./

 data (vwave(i),i=193,224)/8075.,8100.,8125.,8150.,8175.,8200.,8225.,8250.,
             8275.,8300.,8325.,8350.,8375.,8400.,8425.,8450.,
             8475.,8500.,8525.,8550.,8575.,8600.,8625.,8650.,
             8675.,8700.,8725.,8750.,8775.,8800.,8825.,8850./

 data (vwave(i),i=225,256)/8875.,8900.,8925.,8950.,8975.,9000.,9025.,9050.,
             9075.,9100.,9125.,9150.,9175.,9200.,9225.,9250.,
             9275.,9300.,9325.,9350.,9375.,9400.,9425.,9450.,
             9475.,9500.,9525.,9550.,9575.,9600.,9625.,9650./

 data (vwave(i),i=257,288)/9675.,9700.,9725.,9750.,9775.,9800.,9825.,9850.,
             9875.,9900.,9925.,9950.,9975.,10000.,10025.,10050.,
             10075.,10100.,10125.,10150.,10175.,10200.,10225.,10250.,
             10275.,10300.,10325.,10350.,10375.,10400.,10425.,10450./

 data (vwave(i),i=289,291)/10475.,10500.,12000./


 data (vega(i),i=1,32)/0.089,0.089,0.109,0.122,0.124,0.126,0.137,0.15,
             0.161,0.17,0.173,0.176,0.184,0.189,0.194,0.19,
             0.186,0.183,0.145,-0.004,-0.164,-0.329,-0.46,-0.611,
             -0.723,-0.707,-0.824,-0.738,-0.684,-0.864,-0.949,-0.938/

 data (vega(i),i=33,64)/-0.771,-0.633,-0.729,-0.85,-0.837,-0.822,-0.81,-0.789,
             -0.756,-0.662,-0.474,-0.477,-0.601,-0.673,-0.66,-0.644,
             -0.626,-0.609,-0.585,-0.541,-0.549,-0.539,-0.525,-0.509,
             -0.488,-0.468,-0.448,-0.426,-0.409,-0.388,-0.303,-0.135/

 data (vega(i),i=65,96)/-0.14,-0.258,-0.318,-0.304,-0.284,-0.269,-0.252,-0.235,
             -0.221,-0.205,-0.19,-0.174,-0.159,-0.143,-0.13,-0.115,
             -0.102,-0.086,-0.077,-0.063,-0.048,-0.035,-0.02,-0.005,
             0.009,0.027,0.042,0.055,0.064,0.08,0.099,0.115/

 data (vega(i),i=97,128)/0.13,0.146,0.157,0.172,0.187,0.202,0.21,0.229,
             0.244,0.26,0.276,0.292,0.306,0.323,0.337,0.35,
             0.361,0.373,0.388,0.403,0.418,0.433,0.448,0.464,
             0.473,0.485,0.5,0.51,0.516,0.524,0.541,0.555/

 data (vega(i),i=129,160)/0.57,0.59,0.688,0.795,0.84,0.723,0.639,0.651,
             0.663,0.675,0.684,0.696,0.709,0.722,0.737,0.753,
             0.765,0.78,0.794,0.8,0.812,0.819,0.832,0.846,
             0.86,0.872,0.887,0.902,0.912,0.922,0.932,0.943/

 data (vega(i),i=161,192)/0.957,0.972,0.983,0.995,1.004,1.012,1.024,1.036,
             1.047,1.058,1.071,1.08,1.095,1.103,1.118,1.126,
             1.135,1.152,1.16,1.168,1.184,1.193,1.201,1.209,
             1.218,1.227,1.245,1.254,1.264,1.273,1.282,1.291/

 data (vega(i),i=193,224)/1.301,1.311,1.32,1.33,1.351,1.361,1.371,1.381,
             1.393,1.404,1.414,1.426,1.436,1.448,1.47,1.481,
             1.493,1.505,1.529,1.541,1.566,1.592,1.541,1.482,
             1.436,1.437,1.591,1.63,1.493,1.448,1.566,1.671/

 data (vega(i),i=225,256)/1.541,1.493,1.493,1.505,1.616,1.833,1.833,1.656,
             1.59,1.603,1.629,1.726,1.833,1.951,1.969,1.849,
             1.684,1.643,1.656,1.67,1.684,1.684,1.712,1.755,
             1.848,1.987,2.102,2.103,1.915,1.801,1.755,1.741/

 data (vega(i),i=257,288)/1.755,1.755,1.786,1.816,1.832,1.847,1.864,1.863,
             1.864,1.879,1.896,1.948,2.003,2.08,2.208,2.327,
             2.12,2.06,2.002,1.982,1.964,1.946,1.964,1.964,
             1.98,1.979,1.998,1.997,2.014,2.013,2.032,2.031/

 data (vega(i),i=289,291)/2.049,2.049,2.049/

 begin
	#INITIALIZE local variables
	status	=	true
	H = 6.62620E-27
	C = 2.997925E+18
	call strcpy ("PHOTLAM", forms[1,1], 7)
	call strcpy ("photlam", forms[1,2], 7)
	call strcpy ("COUNTS", forms[1,3], 7)
	call strcpy ("counts", forms[1,4], 7)
	call strcpy ("flam", forms[1,5], 4)
	call strcpy ("FLAM", forms[1,6], 4)
	call strcpy ("fnu", forms[1,7], 3)
	call strcpy ("FNU", forms[1,8], 3)
	call strcpy ("photnu", forms[1,9], 6)
	call strcpy ("PHOTNU", forms[1,10], 6)
	call strcpy ("jy", forms[1,11], 2)
	call strcpy ("JY", forms[1,12], 2)
	call strcpy ("mjy", forms[1,13], 3)
	call strcpy ("MJY", forms[1,14], 3)
 	call strcpy ("abmag", forms[1,15], 5)
	call strcpy ("ABMAG", forms[1,16], 5)
	call strcpy ("stmag", forms[1,17], 5)
 	call strcpy ("STMAG", forms[1,18], 5)
	call strcpy ("vegamag", forms[1,19], 7)
	call strcpy ("VEGAMAG", forms[1,20], 7)
	call strcpy ("obmag", forms[1,21], 5)
	call strcpy ("OBMAG", forms[1,22], 5)

	if ( nwave > 0 ){
           for ( i = 1; i <= TWONFORM; i=i+1 ){
              if ( streq ( inform, forms[1,i] )){
                 inform_flag = i
              }
              if ( streq (outform, forms[1,i] )){
                 outform_flag = i
              }
           }
	   if ( streq(inform, outform )){
	      do i = 1, nwave {
	         specout[i] = specin[i]
	      }
	   }else{
	      switch (inform_flag) {
	         case 1,2:
	            do i = 1, nwave {
		       specout[i] = specin[i]
	            }
		 case 3,4:
		    do i = 1, nwave {
		       if ( !IS_INDEFR (specin[i]) ){ 
		          im = MAX(i-1,1)
		          ip = MIN(i+1,nwave)
		          specout[i] = specin[i] / 
		               ( ABS( wave[ip]-wave[im] ) * 0.5 )
		       }else{
                          specout[i] = specin[i]
		       }
		    }
		 case 5,6:
		    factor = 1. / ( H * C )
		    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
		          specout[i] = factor * specin[i] * wave[i]
                       }else{
			  specout[i] = specin[i]
                       }
		    }
                 case 7,8:
                    factor = 1. / H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
                          specout[i] = factor * specin[i] / wave[i]
                       }else{
                          specout[i] = specin[i]
                       }
                    }
                 case 9,10:
                    factor = C
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
                          specout[i] = factor * specin[i] / 
                                         ( wave[i] * wave[i] )
                       }else{
                          specout[i] = specin[i]
                       }
                    }
                 case 11,12:
                    factor = 1.E-23 / H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
                          specout[i] = factor * specin[i] / wave[i]
                       }else{
                          specout[i] = specin[i]
                       }
                    }
                 case 13,14:
                    factor = 1.E-26 / H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
                          specout[i] = factor * specin[i] / wave[i]
                       }else{
                          specout[i] = specin[i]
                       }
                    }
		 case 15,16:
                    factor = 1. / H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
                          specout[i] = factor / wave[i] *
                               10.**( -0.4 * ( specin[i] - ABZERO ))
                       }else{
                          specout[i] = specin[i]
                       }
                    }
                 case 17,18:
                    factor = 1. / ( H * C )
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specin[i]) ){
                          specout[i] = factor * wave[i] *
                               10.**( -0.4 * (specin[i] - STZERO ))
                       }else{
                          specout[i] = specin[i]
                       }
                    }
		 case 19,20:
		    call malloc ( vdata, nwave, TY_REAL )
		    call linterp ( nvega, vwave, vega, 
					nwave, wave, Memr(vdata) )
                    factor = 1. / ( H * C )
                    do i = 1, nwave {
                       if ( Memr[vdata+i-1] != INDEFR ){
                          Memr[vdata+i-1] = factor * wave[i] *
                               10.**( -0.4 * (Memr[vdata+i-1] - STZERO ))
		          if ( specin[i] == 0.0 ){
                             specout[i] = Memr[vdata+i-1]
		          }else if ( IS_INDEFR (specin[i]) ){
			     specout[i] = INDEFR
		          }else{
			     specout[i] = Memr[vdata+i-1] *
		 		   10.0 ** ( -0.4 * specin[i] )
			  }
                       }else{
		          specout[i] = INDEFR
                       }
                    }
		    call mfree ( vdata, TY_REAL )
	         case 21,22:
	            do i = 1, nwave {
		       if ( !IS_INDEFR (specin[i]) ){ 
		          im = MAX(i-1,1)
		          ip = MIN(i+1,nwave)
                          specout[i] = 10.**(-0.4 * specin[i]) / 
	                               ( ABS( wave[ip]-wave[im] ) * 0.5 )
		       }else
                          specout[i] = specin[i]
		    }

	         default:
                    status = false
              }
	      switch (outform_flag) {
	         case 1,2:
                    return
		 case 3,4:
		    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
		          im = MAX(i-1,1)
		          ip = MIN(i+1,nwave)
		          specout[i] = specout[i] *
		               ( ABS( wave[ip]-wave[im] ) * 0.5 )
                       }
		    }
		 case 5,6:
		    factor = H * C
		    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
		          specout[i] = factor * specout[i] / wave[i]
                       }
		    }
                 case 7,8:
                    factor = H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
                          specout[i] = factor * specout[i] * wave[i]
                       }
                    }
                 case 9,10:
                    factor = 1. / C
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
                          specout[i] = factor * specout[i] * 
                                          wave[i] * wave[i]
                       }
                    }
                 case 11,12:
                    factor = 1.E+23 * H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
                          specout[i] = factor * specout[i] * wave[i]
                       }
                    }
                 case 13,14:
                    factor = 1.E+26 * H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
                          specout[i] = factor * specout[i] * wave[i]
                       }
                    }
		 case 15,16:
                    factor = H
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
                          specout[i] = factor * specout[i] * wave[i]
                          if ( specout[i] > 0. ){
                             specout[i] = -2.5 * ALOG10( specout[i] ) + 
                                                   ABZERO
                          }else{
                             specout[i] = 1000.
                          }
                       }
                    }
                 case 17,18:
                    factor = H * C
                    do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
                          specout[i] = factor * specout[i] / wave[i]
                          if ( specout[i] > 0. ){
                             specout[i] = -2.5 * ALOG10( specout[i] ) + 
							STZERO
                          }else{
                             specout[i] = 1000.
                          }
                       }
                    }
		 case 19,20:
		    call malloc ( vdata, nwave, TY_REAL )
		    call linterp ( nvega, vwave, vega, 
					nwave, wave, Memr(vdata) )
                    factor = 1. / ( H * C )
                    do i = 1, nwave {
                       if ( Memr[vdata+i-1] != INDEFR ){
                          Memr[vdata+i-1] = factor * wave[i] *
                               10.**( -0.4 * (Memr[vdata+i-1] - STZERO ))
			  if ( !IS_INDEFR (specout[i]) ) {
 		             if ( Memr[vdata+i-1] > 0.0 && specout[i] > 0.0){
                                specout[i] = 
			         -2.5 * ALOG10( specout[i] / Memr[vdata+i-1] )
		             }else{
		      	        specout[i] = 100.0
                             }
                          }
		       }else{
			  specout[i] = INDEFR
                       }
                    }
		    call mfree ( vdata, TY_REAL )
	         case 21,22:
	            do i = 1, nwave {
                       if ( !IS_INDEFR (specout[i]) ){
		          im = MAX(i-1,1)
		          ip = MIN(i+1,nwave)
	                  if ( specout[i] > 0. ) {
		             specout[i] = -2.5 * alog10( specout[i] *
		                        ( ABS( wave[ip]-wave[im] ) * 0.5 ) )
	                  } else
	                     specout[i] = INDEFR
                       }
		    }

                 default:
                    status = false
              }
           }
           return
	}  
	status = false
        return
end
