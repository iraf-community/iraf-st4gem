# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	<mach.h>

#--------------------------------------------------------------------19 Feb 97--
.help redlaw.x Mar95 nebular/fivel
.ih
NAME
gal_redlaw -- Calculates the ISEF for the Galaxy by Savage & Mathis (1979)
ccm_redlaw -- Calculates the ISEF for the Galaxy by Cardelli, Clayton & Mathis (1989)
lmc_redlaw -- Calculates the ISEF for the LMC by Howarth (1983)  
smc_redlaw -- Calculates the ISEF for the SMC by Prevot, et al. (1984)  
jbk_redlaw -- Returns the ISEF of Whitford (1958), as adapted by Kaler (1976)
.ih
DESCRIPTION
The following procedures return an evaluation of the interstellar 
extinction function (ISEF) as published by various authors.  Often 
the ISEF is expressed as A(lambda)/E(B-V), and is normalized such 
that the function is zero at the V passband, and 1.00 at B.  These 
functions are expressed in the literature as tabulated functions, 
polynomial expressions, or both.  The parameter is wavelength in 
inverse microns, and the result is color excess in magnitudes.  The 
functions are defined or extended to about lambda=1000 Ang; the form 
of the ISEF is not known shortward of that.  

The routines here return A(lambda)/A(V) by renormaling the extinction 
function to 0.0 at zero inverse wavelength, and 1.0 at V.  The ratio 
of total to selective absorption in the V passband, R_v, is thought to 
average 3.10 for most lines of sight in the Galaxy and in the LMC; it 
may be somewhat lower in the SMC.  It is the responsibility of the 
calling function to re-normalize (i.e., multiply) the returned value 
by the choice of R_v.  Note that R_v is not a parameter for any of 
the ISEF's here except for CCM.  

.endhelp
#-------------------------------------------------------------------------------
#  Macro for linear interpolation/extrapolaton:
define	LIN_INTERP ( (($2)-($1)) / (($4)-($3)) * (($5)-($3)) + ($1))

#-------------------------------------------------------------------------------
# GAL_REDLAW --	Evaluation of the average interstellar extinction function for 
#		the Galaxy published by Savage & Mathis 1979, ARA&A, vol. 17, 
#		73-111.  The published extinction function is tabulated as 
#		A(lambda)/E(B-V); returned value is A(lambda)/A(H-beta).  
#
# 		Extinction:
#		lambda < 0.1 um		same formula as lambda = 0.1.
#		0.1 < lambda < 3.4	linear interpolation of tabulated funct.
#		3.4 < lambda		extrapolate linearly in 1/lam 
#
#  Revision History:
#	13-Mar-95 by RAShaw	Based upon synphot library routine EBMVLFUNC

procedure gal_redlaw (wave, extl, npts)

#  Calling arguments:
real	wave[ARB]	# I: wavelength of interest
real	extl[ARB]	# O: extinction evaluation array
int	npts		# I: size of evaluation/wave arrays

#  Local variables:
real 	extab[NTAB]	# tabulated extinction at E(B-V)=1.
int	i, pix		# loop indexes
real	val		# intermediate value
double	x 		# 1. / (wave, in microns)
double	xtab[NTAB]	# tabulated inverse wavelengths

define	NTAB		28	#
define	MIN_WAVE	1000.	# Minimum wave where function is defined.

#  Tabulated inverse wavelengths in microns:
data 	(xtab(i),i=1,NTAB)  / 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82, 
	 2.27,  2.50,  2.91,  3.65,  4.00,  4.17,  4.35,  4.57,  4.76,  5.00, 
	 5.26,  5.56,  5.88,  6.25,  6.71,  7.18,  8.00,  8.50,  9.00,  9.50, 
	10.00 /

#  Tabulated extinction function, A(lambda)/E(B-V):
data	(extab(i),i=1,NTAB) / 0.00,  0.16,  0.38,  0.87,  1.50,  2.32,  3.10, 
	 4.10,  4.40,  4.90,  6.20,  7.29,  8.00,  8.87,  9.67,  9.33,  8.62, 
	 8.00,  7.75,  7.87,  8.12,  8.15,  8.49,  9.65, 10.55, 11.55, 12.90, 
	14.40 /

begin
	do pix = 1, npts {
	    if (wave[pix] < MIN_WAVE || IS_INDEFR(wave[pix]) )
		call error (1, "GAL_REDLAW: Invalid wavelength")

	    # Convert wavelength in angstroms to 1/microns
	    x = 10000.D+0 / wave[pix]

	    # Linear interpolation of extinction law in 1/lam 
	    for (i=2; i <= NTAB && x > xtab[i]; i=i+1) 
	    	;
	    val = LIN_INTERP (extab[i-1], extab[i], xtab[i-1], xtab[i], x) 

	    # Renormalize extinction function to A(lambda)/A(V)
	    extl[pix] = val / 3.1
	}
end


#-------------------------------------------------------------------------------
#  CCM_REDLAW - Computes the interstellar extinction function A(lambda)/A(V) 
#		of Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256. 	
#
#	18 May 93 by RAShaw:	Initial implementation, based upon CCM module 
#				in onedspec.deredden

procedure ccm_redlaw (wave, extl, npts, rv)

#  Calling argument:
real	wave[ARB]	# I: wavelength of emission line, Angstroms
real	extl[ARB]	# O: extinction evaluation array
int	npts		# I: size of evaluation/wave arrays
real	rv		# I: ratio of total to selective extinction

#  Declarations:
real	a, b		# coefficients of extinction function
int	pix		# loop index
double	x		# wavelength in inverse microns
double	y		# generic value

define	MIN_WAVE	1000.	# Minimum wave where function is defined.
begin
	do pix = 1, npts {
	    if (wave[pix] < MIN_WAVE || IS_INDEFR(wave[pix]) )
		call error (1, "CCM_REDLAW: Invalid wavelength")

	    # Convert input wavelength to inverse microns
	    x = 10000.D+0 / wave[pix]

	    # For wavelengths longward of the L passband, linearly interpolate 
	    # to 0. at 1/lambda = 0.  (a=0.08, b=-0.0734 @ x=0.29)
	    if (x < 0.29D+0) {
		a =  0.2759 * x
		b = -0.2531 * x

	    # Compute a(x) and b(x)
	    } else if (x < 1.1D+0) {
	    	y = x ** 1.61
	    	a =  0.574 * y
	    	b = -0.527 * y

	    } else if (x < 3.3D+0) {
	    	y = x - 1.82
	    	a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 +
		    y * (0.72085 + y * (0.01979 + y * (-0.77530 + 
		    y * 0.32999))))))
	    	b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * 
		    (-5.38434 + y * (-0.62251 + y * (5.30260 - y * 2.09002))))))

	    } else if (x < 5.9D+0) {
	    	y = (x - 4.67) ** 2
	    	a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
	    	b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)

	    } else if (x < 8.0D+0) {
	    	y = (x - 4.67) ** 2
	    	a = 1.752 - 0.316 * x - 0.104 / (y + 0.341)
	    	b = -3.090 + 1.825 * x + 1.206 / (y + 0.263)

	    	y = x - 5.9
	    	a = a - 0.04473 * y**2 - 0.009779 * y**3
	    	b = b + 0.2130 * y**2 + 0.1207 * y**3

	    # Truncate range of ISEF to that for 1000 Ang.
	    } else if (x <= 10.0D+0) {
	    	x = min (x, 10.0D+0)
	    	y = x - 8.D+0
	    	a = -1.072 - 0.628 * y + 0.137 * y**2 - 0.070 * y**3
	    	b = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374 * y**3
	    }

	    # Compute A(lambda)/A(V)
	    y = a + b / rv

	    extl[pix] = y
	}
end


#-------------------------------------------------------------------------------
#  LMC_REDLAW -	Evaluation of LMC extinction curve as published by Howarth 
#		1983, MNRAS, 203, 301.  Functional fits from paper are 
#		A(lambda)/E(B-V); returned value is A(lambda)/A(V). 
#
#  Revision History:
#	18-Oct-94 by RAShaw	Initial implementation
#	14-Mar-95 by RAShaw	Return A(lambda)/A(V) instead

procedure lmc_redlaw (wave, extl, npts)

#  Calling arguments:
real	wave[ARB]	#I: evaluation wavelength array
real	extl[ARB]	#O: extinction evaluation array
int	npts		#I: size of evaluation/wave arrays

#  Local variables:
double	delt		# work variable
real 	extab[NTAB]	# tabulated extinction at E(B-V)=1.
int	i, pix		# loop index
real	val		# intermediate value
real	x 		# 1. / (wave, in microns)
real	xtab[NTAB]	# tabulated inverse wavelengths

define	NTAB		7
define	MIN_WAVE	1000.		# Minimum wave where function is defined.

# Tabulated inverse wavelengths in microns:
data 	(xtab(i),i=1,NTAB)  / 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82 /

# Tabulated extinction function, A(lambda)/E(B-V), from Savage & Mathis:
data	(extab(i),i=1,NTAB) / 0.00,  0.16,  0.38,  0.87,  1.50,  2.32,  3.10 /

begin
	do pix = 1, npts {
	    if (wave[pix] < MIN_WAVE || IS_INDEFR(wave[pix]) )
		call error (1, "LMC_REDLAW: Invalid wavelength")

	    # Convert wavelength in angstroms to 1/microns
	    x = 10000. / wave[pix]

	    # Infrared - optical: linear interpolation of Savage & Mathis 1979
	    if ( x <= 1.82) {
	    	for (i=2; i <= NTAB && x > xtab[i]; i=i+1) 
	    	    ;
	    	val = LIN_INTERP (extab[i-1], extab[i], xtab[i-1], xtab[i], x) 

	    # The following polynomial evaluations assume R = 3.1
	    #   Violet
	    } else if ( x <= 2.75 ) { 
	    	delt = x - 1.82
	    	val = 3.1 + (2.04 + 0.094 * delt) * delt
		
	    #   Ultraviolet out to lambda = 1000 A
	    } else { 
	    	x    = min (x, 10.0)
	    	delt = x - 4.557
	    	val  = 3.1 - 0.236 + 0.462 * x + 0.105 * x * x +
		       0.454 / (delt**2 + 0.293)
	    }

	    # Renormalize extinction function to A(lambda)/A(V)
	    extl[pix] = val / 3.1
	}
end


#-------------------------------------------------------------------------------
#  SMC_REDLAW -	Computes the interstellar extinction function of Prevot et al. 
#		(1984), A&A, 132, 389-392.  Tabulated values from paper are 
#		E(lambda)/E(B-V); returned value is A(lambda)/A(v).  
#
#	20-Sep-94 by RAShaw:	Initial implementation
#	14-Mar-95 by RAShaw	Return A(lambda)/A(V) instead

procedure smc_redlaw (wave, extl, npts)

#  Calling argument:
real	wave[ARB]	# I: wavelength of emission line, Angstroms
real	extl[ARB]	# O: extinction evaluation array
int	npts		# I: size of evaluation/wave arrays

#  Local variables:
real	extab[NTAB]	# tabulated extinction function from Prevot
int	i, pix		# loop indexes
real	val		# generic
real	x		# wavelength in inverse microns
real	xtab[NTAB]	# tabulated inverse wavelengths from Prevot

define	NTAB		33
define	MIN_WAVE	1000.		# Minimum wave where function is defined.

# Tabulated inverse wavelengths in microns:
data 	(xtab(i),i=1,NTAB)  / 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82, 
 	 2.35,  2.70,  3.22,  3.34,  3.46,  3.60,  3.75,  3.92,  4.09,  4.28, 
	 4.50,  4.73,  5.00,  5.24,  5.38,  5.52,  5.70,  5.88,  6.07,  6.27, 
	 6.48,  6.72,  6.98,  7.23,  7.52,  7.84 /

# Tabulated extinction function, E(lambda-V)/E(B-V):
data	(extab(i),i=1,NTAB) /-3.10, -2.94, -2.72, -2.23, -1.60, -0.78,  0.00, 
	 1.00,  1.67,  2.29,  2.65,  3.00,  3.15,  3.49,  3.91,  4.24,  4.53, 
     	 5.30,  5.85,  6.38,  6.76,  6.90,  7.17,  7.71,  8.01,  8.49,  9.06,  
	 9.28,  9.84, 10.80, 11.51, 12.52, 13.54 /

begin
	do pix = 1, npts {
	    if (wave[pix] < MIN_WAVE || IS_INDEFR(wave[pix]) )
		call error (1, "SMC_REDLAW: Invalid wavelength")

	    # Convert wavelength in angstroms to 1/microns
	    x = 10000. / wave[pix]
	    x = min (x, 10.0)

	    # Linearly interpolate extinction law in 1/lam 
	    for (i=2; i <= NTAB && x > xtab[i]; i=i+1) 
	    	;
	    val = LIN_INTERP (extab[i-1], extab[i], xtab[i-1], xtab[i], x) 

	    # Renormalize extinction function to A(lambda)/A(V)
	    extl[pix] = val / 3.1 + 1.
	}
end


#-------------------------------------------------------------------------------
#  JBK_REDLAW -	Computes the interstellar extinction function of Whitford 
#		(1958), extended to the UV by Seaton (1977), as adapted by 
#		Kaler (1976). 
#
#	13 May 93 by RAShaw:	initial implementation

procedure jbk_redlaw (wave, extl, npts)

#  Calling argument:
real	wave[ARB]	#I: wavelength of emission line, Angstroms
real	extl[ARB]	#O: extinction evaluation array
int	npts		#I: size of evaluation/wave arrays

#  Local variables:
int	i, pix		# loop indexes
real	x		# wavelength in inverse microns
real	x1, x2		# bounding reference wavelengths, in micron^-1
real	extab[NTAB]	# tabulated extinction function from Kaler (1976)
real	refw[NTAB]	# reference wavelengths

define	NTAB		95	# 
define	MIN_WAVE	1000.	# Minimum wave where function is defined.

# Tabulated wavelengths, Angstroms:
data	refw   /1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550.,
	 1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000., 2050., 
	 2100., 2150., 2200., 2250., 2300., 2350., 2400., 2450., 2500., 2550., 
	 2600., 2650., 2700., 2750., 2800., 2850., 2900., 2950., 3000., 3050., 
	 3100., 3333., 3500., 3600., 3700., 3800., 3900., 4000., 4100., 4200., 
	 4300., 4400., 4500., 4600., 4700., 4800., 4861.3, 5000., 5100., 5200., 
	 5300., 5400., 5500., 5600., 5700., 5800., 5900., 6000., 6100., 6200., 
	 6300., 6400., 6500., 6600., 6700., 6800., 6900., 7000., 7200., 7400., 
	 7600., 7800., 8000., 8200., 8400., 8600., 8800., 9000., 9500., 10000.,
	11000.,12000.,14000.,16000.,20000.,1.D+6/

# Tabulated extinction function:
data	extab  /1.96,  1.78,  1.61,  1.49,  1.37,  1.29,  1.24,  1.20,  1.20,
	 1.20,  1.17,  1.13,  1.11,  1.10,  1.12,  1.17,  1.25,  1.35,  1.45,
	 1.53,  1.60,  1.62,  1.52,  1.40,  1.28,  1.17,  1.06,  0.98,  0.9,
	 0.84,  0.77,  0.72,  0.68,  0.64,  0.60,  0.57,  0.53,  0.51,  0.48,
	 0.46,  0.385, 0.358, 0.33,  0.306, 0.278, 0.248, 0.220, 0.195, 0.168,
	 0.143, 0.118, 0.095, 0.065, 0.040, 0.015, 0.000,-0.030,-0.055,-0.078,
	-0.10, -0.121,-0.142,-0.164,-0.182,-0.201,-0.220,-0.238,-0.254,-0.273,
	-0.291,-0.306,-0.321,-0.337,-0.351,-0.365,-0.377,-0.391,-0.416,-0.441,
	-0.465,-0.490,-0.510,-0.529,-0.548,-0.566,-0.582,-0.597,-0.633,-0.663,
	-0.718,-0.763,-0.840,-0.890,-0.960,-1.000/

begin
	do pix = 1, npts {
	    if (wave[pix] < MIN_WAVE || IS_INDEFR(wave[pix]) )
		call error (1, "SMC_REDLAW: Invalid wavelength")

	    for (i=2; i <= NTAB && refw[i] < wave[pix] ; i=i+1) 
	    	;

	    if ( abs (wave[pix] - refw[i]) < EPSILONR )
	    	extl[pix] = extab[i]

	    else {
	    	x  = 10000. / wave[pix]
#	    	x  = min (x, 10.0)
	    	x1 = 10000. / refw[i-1]
	    	x2 = 10000. / refw[i]
	    	extl[pix] = LIN_INTERP (extab[i-1], extab[i], x1, x2, x)
	    }
	}
end


