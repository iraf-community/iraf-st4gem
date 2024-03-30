# DPHOT.H -- Photometry descriptor structure

define	LEN_PHTSTRUCT	8

define	PHT_NWAVE	Memi[$1]	# number of wavelengths in each array
define	PHT_NBAND	Memi[$1+1]	# number of band arrays
define	PHT_WAVE	Memi[$1+2]	# wavelength array
define	PHT_BAND	Memi[$1+3]	# band arrays
define	PHT_NSTAT	Memi[$1+4]	# number of bandpass statistics
define	PHT_STIM	Memi[$1+5]	# effective stimulus array
define	PHT_PIVOT	Memi[$1+6]	# pivot wavelength array
define	PHT_FWHM	Memi[$1+7]	# equivalent fwhm array
