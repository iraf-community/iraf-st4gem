# Column names, units and formats.
#
# Modified 2/22/96 (I Busko) to implement chi-squared minimization.
# Modified 9/14/98 (IB) to fix problem with FITS output tables.
# Modified 8/30/99 (IB) to fix Y2K problem.

define	DUMMY		"foo"		# dummy entry

define	MIN_COLS	10		# minimum number of columns 
					# in output SDAS table

define	DB_CFILE	"file"		# original file name
define	DB_UFILE	""
define	DB_SFILE	25		# size
define	DB_FFILE	""

define	DB_CTIME	"time"		# time of creation
define	DB_UTIME	""
define	DB_STIME	SZ_TIME		# size
define	DB_FTIME	""

define	DB_CFUNC	"function"	# function type
define	DB_UFUNC	""
define	DB_SFUNC	10		# size
define	DB_FFUNC	""

define	DB_CUNIT	"unit"		# physical units
define	DB_UUNIT	""
define	DB_SUNIT	10		# size
define	DB_FUNIT	""

define	DB_CDEGR	"ncoeff"	# number of function coefficients
define	DB_UDEGR	""
define	DB_FDEGR	""

define	DB_CNPTS	"npoints"	# no. of points in fit
define	DB_UNPTS	""
define	DB_FNPTS	""

define	DB_CXMIN	"xmin"		# xmin for linear function
define	DB_UXMIN	""		# normalization
define	DB_FXMIN	""

define	DB_CXMAX	"xmax"		# xmax for linear function
define	DB_UXMAX	""		# normalization
define	DB_FXMAX	""

define	DB_CCHI		"chisq"		# reduced chi-square of fit
define	DB_UCHI		""
define	DB_FCHI		""

define	DB_CRMS		"rms"		# rms of fit
define	DB_URMS		""
define	DB_FRMS		""

define	DB_CCOEF	"coeff%d"	# function coefficient
define	DB_UCOEF	""
define	DB_FCOEF	""

define	DB_CERR		"err%d"		# coefficient error
define	DB_UERR		""
define	DB_FERR		""


