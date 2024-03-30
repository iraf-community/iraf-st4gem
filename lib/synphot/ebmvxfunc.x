define	GAL1		1
define	GAL2		2
define	GAL3		3
define	SMC		4
define  LMC		5
define	XGAL		6

# EBMVXFUNC -- Extended reddening law extinction function

procedure ebmvxfunc (extval, redlaw, nwave, wave, band)

real	extval		# i: extinction value
char	redlaw[ARB]	# i: extinction function name
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
int	ftype, iwave
pointer	sp, func

string	funclist    "gal1,gal2,gal3,smc,lmc,xgal"
string	badredlaw   "Unrecognized redenning law name"

int	word_match()

begin
	# Convert function name to lower case

	call smark(sp)
	call salloc (func, SZ_FNAME, TY_CHAR)

	call strcpy (redlaw, Memc[func], SZ_FNAME)
	call strfix (Memc[func])

	# Call selected extinction function to compute 
	# extinction in magnitudes

	ftype = word_match (Memc[func], funclist)

	switch (ftype) {
	case GAL1:
	    call ebmv_gal1 (nwave, wave, band)
	case GAL2:
	    call ebmv_gal2 (nwave, wave, band)
	case GAL3:
	    call ebmv_gal3 (nwave, wave, band)
	case SMC:
	    call ebmv_smc (nwave, wave, band)
	case LMC:
	    call ebmv_lmc (nwave, wave, band)
	case XGAL:
	    call ebmv_xgal (nwave, wave, band)
	default:
	    call synphoterr (badredlaw, redlaw)
	}

	# Convert extinction from magnitudes

	do iwave = 1, nwave
	    band[iwave] = 10.0 ** (-0.4 * extval * band[iwave])

	call sfree (sp)
end

# EBMV_GAL1 -- First galactic reddening law (Seaton's)
#
# Seaton's paper in M.N.R.A.S. vol 187, page 75p (1979).  
# the formulae are based on an adopted value of R = 3.20.
#
# Note that Seaton's representation of of the interstellar reddening law
# differs substantially from Schild's representation (Astron. J. 82, 339,
# table ii, 1977) in the region of overlap.  Schild also adopted r = 3.20.
# For wavelengths > 3704 angstroms, the program interpolates 
# linearly in 1/lambda in Seaton's table 3. For wavelengths < 3704 
# angstroms, the program uses the formulae from Seaton's table 2. 
# The formulae match at the endpoints of their respective intervals. 
# there is a mismatch of 0.009 mag/ebmv at nu=2.7 (lambda=3704 angstroms).
# Seaton's tabulated value of 1.44 mags at 1/lambda = 1.1 may be in error;
# 1.64 seems more consistent with his other values. 
#
# Sources:
#	lambda < 1000		same formula as lambda = 1000.
#	1000 < lambda < 3704	Seaton(1979) MNRAS 187,73p.
#	3704 < lambda < 10,000	Nandy(1975) A+A 44, 195. (corrected to R=3.2)
#	10000 < lambda		extrapolate linearly in 1/lam (can be improved)

procedure ebmv_gal1 (nwave, wave, band)

int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
real	xtable[19]	# tabulated inverse wavelengths
real 	etable[19]	# tabulated extinctin at E(B-V)=1.

int	i, iwave
real	x

data 	(xtable(i),i=1,19) / 0.,  1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
			     1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
			     2.2, 2.3, 2.4, 2.5, 2.6, 2.7 /

data	(etable(i),i=1,19) / 0.,   1.36, 1.64, 1.84, 2.04, 2.24, 2.44,
                             2.66, 2.88, 3.14, 3.36, 3.56, 3.77,
     			     3.96, 4.15, 4.26, 4.40, 4.52, 4.64 /

begin

	do iwave = 1, nwave {
	    # Convert wavelength in angstroms to 1/microns
	    x = 10000.0 / wave[iwave]

	    # Infrared: extend optical results linearly to 0 at 1/lam = 0
	    if (x <= 1.0) {
		band[iwave] = etable[2] * x * x

	    # Optical region interpolates in Seaton's table 3
	    } else if (x < 2.7) {
		call ebmv_interp (19, xtable, etable, x, band[iwave])

	    # Ultraviolet uses analytic formulae from Seaton's table 2
	    } else if (x < 3.65) {
		band[iwave] = 1.56 + 1.048 * x + 1.01 /
			 ((x - 4.6) * (x - 4.6) + 0.280)

	    # More ultraviolet
	    } else if (x < 7.14) {
		band[iwave] = 2.29 + 0.848 * x + 1.01 /
			 ((x - 4.6) * (x - 4.6) + 0.280)

	    # And more ultraviolet still
	    } else {
		x = min(x, 50.0)
		band[iwave] = 16.17 + x * (-3.20 + 0.2975 * x)
	    }
	}

end

# EBMV_GAL2 -- Second galactic reddening law (Savage and Mathis's)
#
# Evaluation of the average interstellar extinction function the Galaxy 
# published by Savage & Mathis 1979, ARA&A, vol. 17, 73-111.

procedure ebmv_gal2 (nwave, wave, band)

int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
real	xtable[28]	# tabulated inverse wavelengths
real 	etable[28]	# tabulated extinction at E(B-V)=1.

int	i, iwave
real	x

#  Tabulated inverse wavelengths in microns:
data 	(xtable(i),i=1,28)  / 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82, 
	 2.27,  2.50,  2.91,  3.65,  4.00,  4.17,  4.35,  4.57,  4.76,  5.00, 
	 5.26,  5.56,  5.88,  6.25,  6.71,  7.18,  8.00,  8.50,  9.00,  9.50, 
	10.00 /

#  Tabulated extinction function, A(lambda)/E(B-V):
data	(etable(i),i=1,28)  / 0.00,  0.16,  0.38,  0.87,  1.50,  2.32,  3.10, 
	 4.10,  4.40,  4.90,  6.20,  7.29,  8.00,  8.87,  9.67,  9.33,  8.62, 
	 8.00,  7.75,  7.87,  8.12,  8.15,  8.49,  9.65, 10.55, 11.55, 12.90, 
	14.40 /

begin
	do iwave = 1, nwave {
	    # Convert wavelength in angstroms to 1/microns
	    x = 10000.0 / wave[iwave]

	    # Interpolate for extinction
	    call ebmv_interp (28, xtable, etable, x, band[iwave])
	}
end

# EBMV_GAL3 -- Third galactic reddening law (Cardelli, Clayton & Mathis)
#
# Computes the interstellar extinction function A(lambda)/A(V) 
# of Cardelli, Clayton & Mathis 1989, ApJ 345, 245-256. 	

procedure ebmv_gal3 (nwave, wave, band)

int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
int	iwave
real	x, y, a, b

begin
	do iwave = 1, nwave {
	    # Convert input wavelength to inverse microns
	    x = 10000.0 / wave[iwave]

	    # Linear continuation of formula for 0.3 <= x < 1.1
	    if (x < 0.30) {
		a =  0.2754 * x
		b = -0.2528 * x

	    # Compute a(x) and b(x)
	    } else if (x < 1.1) {
	    	y = x ** 1.61
	    	a =  0.574 * y
	    	b = -0.527 * y

	    } else if (x < 3.3) {
	    	y = x - 1.82
	    	a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 +
		    y * (0.72085 + y * (0.01979 + y * (-0.77530 + 
		    y * 0.32999))))))
	    	b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * 
		    (-5.38434 + y * (-0.62251 + y * (5.30260 - y * 2.09002))))))

	    } else if (x < 5.9) {
	    	a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
	    	b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)

	    } else if (x < 8.0) {
	    	a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
	    	b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)

	    	y = x - 5.9
	    	a = a - 0.04473 * y**2 - 0.009779 * y**3
	    	b = b + 0.2130 * y**2 + 0.1207 * y**3

	    # Truncate range of ISEF to that for 1000 Ang.
	    } else {
	    	x = min (x, 10.0)
	    	y = x - 8.0
	    	a = -1.073 - 0.628 * y + 0.137 * y**2 - 0.070 * y**3
	    	b = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374 * y**3
	    }

	    # Compute A(lambda)
	    band[iwave] = 3.1 * a + b
	}

end

# EBMV_SMC -- Reddening law for Small Magellanic Cloud
#
# Computes the interstellar extinction function of Prevot et al. 
# (1984), A&A, 132, 389-392. 

procedure ebmv_smc (nwave, wave, band)

int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
real	xtable[33]	# tabulated inverse wavelengths
real 	etable[33]	# tabulated extinction at E(B-V)=1.

int	i, iwave
real	x, val

# Tabulated inverse wavelengths in microns:
data 	(xtable(i),i=1,33)  / 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.82, 
 	 2.35,  2.70,  3.22,  3.34,  3.46,  3.60,  3.75,  3.92,  4.09,  4.28, 
	 4.50,  4.73,  5.00,  5.24,  5.38,  5.52,  5.70,  5.88,  6.07,  6.27, 
	 6.48,  6.72,  6.98,  7.23,  7.52,  7.84 /

# Tabulated extinction function, E(lambda-V)/E(B-V):
data	(etable(i),i=1,33)  /-3.10, -2.94, -2.72, -2.23, -1.60, -0.78,  0.00, 
	 1.00,  1.67,  2.29,  2.65,  3.00,  3.15,  3.49,  3.91,  4.24,  4.53, 
     	 5.30,  5.85,  6.38,  6.76,  6.90,  7.17,  7.71,  8.01,  8.49,  9.06,  
	 9.28,  9.84, 10.80, 11.51, 12.52, 13.54 /

begin
	do iwave = 1, nwave {
	    # Convert wavelength in angstroms to 1/microns
	    x = 10000.0 / wave[iwave]
	    x = min (x, 10.0)

	    # Interpolate for extinction
	    call ebmv_interp (33, xtable, etable, x, val)

	    band[iwave] = val + 3.1
	}
end

# EBMV_LMC -- Reddening law for the Large Magellanic Cloud
#
# Evaluation of LMC extinction curve as published by Howarth 1983, 
# MNRAS, 203, 301. 

procedure ebmv_lmc (nwave, wave, band)

int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
real	xtable[7]	# tabulated inverse wavelengths
real 	etable[7]	# tabulated extinction at E(B-V)=1.

int	i, iwave
real	x, delt

# Tabulated inverse wavelengths in microns:
data 	(xtable(i),i=1,7) / 0.00,  0.29,  0.45,  0.80,  1.11,  1.43,  1.83 /

# Tabulated extinction function, A(lambda)/E(B-V), from Savage & Mathis:
data	(etable(i),i=1,7) / 0.00,  0.16,  0.38,  0.87,  1.50,  2.32,  3.10 /

begin
	do iwave = 1, nwave {
	    # Convert wavelength in angstroms to 1/microns
	    x = 10000.0 / wave[iwave]

	    # Interpolate for extinction
	    if ( x < 1.83) {
		call ebmv_interp (7, xtable, etable, x, band[iwave])

	    #   Violet
	    } else if ( x <= 2.75 ) { 
	    	delt = x - 1.83
	    	band[iwave] = 3.1 + (2.04 + 0.094 * delt) * delt
		
	    #   Ultraviolet out to lambda = 1000 A
	    } else { 
	    	x    = min (x, 10.0)
	    	delt = x - 4.557
	    	band[iwave]  = 3.1 - 0.236 + 0.462 * x + 0.105 * x * x +
			       0.454 / (delt**2 + 0.293)
	    }
	}
end

# EBMV_XGAL -- Extragalactic reddening function
#
# Computes the extragalactic extinction function 
# of Calzetti, Kinney and Storchi-Bergmann, 1994, ApJ, 429, 582 

procedure ebmv_xgal (nwave, wave, band)

int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
int	iwave
real	x, y

begin
	do iwave = 1, nwave {
	    # Convert input wavelength to inverse microns
	    x = 10000.0 / wave[iwave]

	    # Formula from paper with zero point moved to (x = 0)
	    y = ((0.011 * x - 0.198) * x + 1.509) * x

	    # Normalize the result according to Kailash Sahu's calculations
	    band[iwave] = 2.43 * y
	}
end

# EBMV_INTERP -- Linear interpolation routine for ebmv functions

procedure ebmv_interp (ntab, xtable, etable, x, ext)

int	ntab		# i: table length
real	xtable[ARB]	# i: inverse wavelength table
real	etable[ARB]	# i: extinction table 
real	x		# i: inverse wavelength (in microns)
real	ext		# o: interpolated extinction value (in magnitudes)
#--
int	itab
real	a, b

begin
	# Starting at 2 and ending at ntab-1 ensures a pair of points
	# for the interpolation even if the input value lies outside
	# of the table and the interpolation is thus an extrapolation

	for (itab = 2; itab < ntab; itab = itab + 1) {
	    if (x <= xtable[itab])
		break
	}

	# Standard linear interpolation formula 

	a = (xtable[itab] - x) / (xtable[itab] - xtable[itab-1])
	b = 1.0 - a

	ext = a * etable[itab-1] + b * etable[itab]
end
