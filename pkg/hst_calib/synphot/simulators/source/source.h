# SOURCE.H -- Spectral source structure

define	SZ_SRCSTRUCT	10
define	SZ_SRCSCALARS	3

define	SRC_NWAVE	Memi[$1]	# length of wavelength array
define	SRC_NXPSF	Memi[$1+1]	# x dimension of PSF
define	SRC_NYPSF	Memi[$1+2]	# y dimension of PSF
define	SRC_NXLSF	Memi[$1+3]	# length of LSF
define	SRC_SCALARS	Memi[$1+4]	# real scalar variables
define	SRC_WAVPTR	Memi[$1+5]	# wavelength array
define	SRC_SPECPTR	Memi[$1+6]	# source spectrum array
define	SRC_PSFPTR	Memi[$1+7]	# psf weight array
define	SRC_LSFPTR	Memi[$1+8]	# lsf weight array
define	SRC_SHAPE	Memi[$1+9]	# source shape

define	SRC_XPOS	Memr[SRC_SCALARS($1)]
define	SRC_YPOS	Memr[SRC_SCALARS($1)+1]
define	SRC_FLUX	Memr[SRC_SCALARS($1)+2]
define	SRC_WAVE	Memr[SRC_WAVPTR($1)]
define	SRC_SPEC	Memr[SRC_SPECPTR($1)]
define	SRC_PSF		Memr[SRC_PSFPTR($1)]
define	SRC_LSF		Memr[SRC_LSFPTR($1)]
