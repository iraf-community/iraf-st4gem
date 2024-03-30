# DSPEC.H -- Spectrum descriptor structure

define	LEN_SPCSTRUCT	7

define	SPC_NWAVE	Memi[$1]	# number of wavelengths in each array
define	SPC_NFLUX	Memi[$1+1]	# number of flux arrays
define	SPC_WAVE	Memi[$1+2]	# wavelength array
define	SPC_FLUX	Memi[$1+3]	# flux arrays
define	SPC_STIM	Memi[$1+4]	# effective stimulus
define	SPC_PIVOT	Memi[$1+5]	# pivot wavelength
define	SPC_FWHM	Memi[$1+6]	# equivalent full width, half max
