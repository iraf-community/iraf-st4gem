# DBAND.H -- Bandpass descriptor structure

define	LEN_BNDSTRUCT	4

define	BND_NWAVE	Memi[$1]	# number of wavelengths in each array
define	BND_NBAND	Memi[$1+1]	# number of bandpass arrays
define	BND_WAVE	Memi[$1+2]	# wavelength array
define	BND_BAND	Memi[$1+3]	# bandpass arrays

