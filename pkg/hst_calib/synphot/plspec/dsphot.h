# DSPHOT.H -- Spectrophotometry descriptor structure

define	LEN_SPTSTRUCT	6

define	SPT_NWAVE	Memi[$1]	# number of wavelengths in each array
define	SPT_NFLUX	Memi[$1+1]	# number of flux arrays
define	SPT_WAVE	Memi[$1+2]	# wavelength array
define	SPT_FLUX	Memi[$1+3]	# flux array
define	SPT_ERROR	Memi[$1+4]	# error array
define	SPT_FWHM	Memi[$1+5]	# equivalent fwhm array


