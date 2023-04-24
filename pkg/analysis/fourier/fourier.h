# FOURIER.H -- Include file for SDAS/IRAF Fourier analysis program.
# Defines two structures:  FT for info about a file; WORK for pointers
# to work space for NCAR and for phase shifting and centering output.

define	FT_MAX_NDIM	2		# Maximum dimension of image
define	MAX_CTPAIRS	30		# Maximum number of CTYPEs

# Description of FT structure.

define	LEN_FT		44		# Length of Fourier structure

define	FT_REAL		Memi[$1]	# Real part exists?
define	FT_IMAG		Memi[$1+1]	# Imaginary part exists?

define	FT_REPT		Memi[$1+2]	# Real image pointer or null
define	FT_IMPT		Memi[$1+3]	# Imaginary image pointer or null
define	FT_IMAGE	Memi[$1+4]	# =FT_REPT, or FT_IMPT if FT_REAL=NO

# These are for the names of the images for the real and imaginary parts.
define	FT_NAME_R_PTR	Memi[$1+5]	# pointer to name of real part
define	FT_NAME_I_PTR	Memi[$1+6]	# pointer to name of imaginary part
define	FT_NAME_R	Memc[FT_NAME_R_PTR($1)] # name of real part
define	FT_NAME_I	Memc[FT_NAME_I_PTR($1)] # name of imaginary part

# Note:  NAXIS is limited to 2; see FT_MAX_NDIM above.
# The macro for the CD matrix is defined as follows:  FT_CD(ft,i,j) = CDi_j

define	FT_NAXIS    Memi[$1+7]			# Number of data axes
#		 same as P2D($1+8) for FT_CRVAL(ft,1)
define	FT_CRVAL    Memd[P2D($1+ 8+($2-1)*2)]	# Value at reference pixel
define	FT_CRPIX    Memd[P2D($1+12+($2-1)*2)]	# Reference pixel
define	FT_CD       Memd[P2D($1+16+($2-1)*2+($3-1)*4)]	# CD matrix
define	FT_CTYPE    Memc[P2C($1+24+($2-1)*10)]	# Coordinate type for axis

# next available location is ($1+44)

# Description of WORK structure.

define	LEN_WORK	9		# length of WORK structure

define	N_FILES		Memi[$1]	# number of input files (1 or 2)
define	NPTS		Memi[$1+1]	# number of points in input array
define	SHIFT		Memi[$1+2]	# amount of shift; used if center=true
define	TRIGTAB		Memi[$1+3]	# pointer to NCAR trig table
define	XWORK		Memi[$1+4]	# pointer to complex array
define	XWORK2		Memi[$1+5]	# pointer; used if N_FILES=2
define	COSTAB		Memi[$1+6]	# pointer to cos table for phase shift
define	SINTAB		Memi[$1+7]	# pointer to sine table
define	C_COPY		Memi[$1+8]	# pointer to space for centering output
