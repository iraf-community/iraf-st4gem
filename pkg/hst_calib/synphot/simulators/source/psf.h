# PSF.H -- PSF descriptor

define LEN_PSFSTRUCT	5

define PSF_NUMBER	Memi[$1]	# Number of psfs
define PSF_NXPIX	Memi[$1+1]	# Number of pixels in each dimension
define PSF_NYPIX	Memi[$1+2]	# Number of pixels in each dimension
define PSF_WAVPTR	Memi[$1+3]	# Pointer to wavelength array
define PSF_BUFFER	Memi[$1+4]	# Psf buffer

define PSF_WAVE		Memr[PSF_WAVPTR($1)+($2)-1]
define PSF_ARYPTR	PSF_BUFFER($1)+(($2)-1)*PSF_NXPIX($1)*PSF_NYPIX($1)
define PSF_ARRAY	Memr[PSF_ARYPTR($1,$2)]

