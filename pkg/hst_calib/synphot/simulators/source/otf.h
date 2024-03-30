# OTF.H -- Optical transfer function (PSF or LSF) descriptor


define MAXOTF		10000	     # maximum size of OTF dimension

define LEN_OTFSTRUCT	5

define OTF_NUMBER	Memi[$1]	# Number of otfs
define OTF_NXPIX	Memi[$1+1]	# Number of pixels in each dimension
define OTF_NYPIX	Memi[$1+2]	# Number of pixels in each dimension
define OTF_WAVPTR	Memi[$1+3]	# Pointer to wavelength array
define OTF_BUFFER	Memi[$1+4]	# OTF buffer

define OTF_WAVE		Memr[OTF_WAVPTR($1)+($2)-1]
define OTF_ARYPTR	OTF_BUFFER($1)+(($2)-1)*OTF_NXPIX($1)*OTF_NYPIX($1)
define OTF_ARRAY	Memr[OTF_ARYPTR($1,$2)]

