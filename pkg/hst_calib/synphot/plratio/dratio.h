# DRATIO.H -- Ratio descriptor structure

define	LEN_RATSTRUCT	4

define	RAT_NWAVE	Memi[$1]	# number of wavelengths in each array
define	RAT_NFLUX	Memi[$1+1]	# number of flux arrays
define	RAT_WAVE	Memi[$1+2]	# wavelength array
define	RAT_FLUX	Memi[$1+3]	# flux arrays

