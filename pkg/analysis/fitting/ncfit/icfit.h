# The ICFIT data structure. Modified to allow layering on top of the
# non-linear 1-d fitting package NLFIT.

include	"../nlfit/nlfit.h"

# The function names must be in the same order as in the definitions 
# in nlfit.h

define	FUNCTIONS "|bbody|powerlaw|composite|twobbody|galprof|galbulge|galdisk|Gaussians|cgauss|user"

# Same as above.

define	METHODS "|amoeba|Levenberg-Marquardt"

#-------------------------------------------------------------------
# Types of y axis

define	AXISTYPE "|linear|log|magnitude"

# Types of y axis
define	IC_LINEAR	0
define	IC_LOG		1
define	IC_MAG		2

# Types of graphic lines.

define	IC_CONTINUOUS	0		# Continuous line
define	IC_BOX		1		# Markers whith error bars

define	IC_NGKEYS	5		# Number of graph keys
define	IC_LENSTRUCT	47		# Length of ICFIT structure

# User fitting parameters
define	IC_FUNCTION	Memi[$1]	# Function type
define	IC_NPAR		Memi[$1+1]	# Number of coefficients
define	IC_SAMPLE	Memi[$1+2]	# Pointer to sample string
define	IC_NAVERAGE	Memi[$1+3]	# Sampling averaging bin
define	IC_NITERATE	Memi[$1+4]	# Number of rejection interation
define	IC_XMIN		Memr[$1+5]	# Minimum value for curve
define	IC_XMAX		Memr[$1+6]	# Maximum value for curve
define	IC_LOW		Memr[$1+7]	# Low rejection value
define	IC_HIGH		Memr[$1+8]	# Low rejection value
define	IC_GROW		Memr[$1+9]	# Rejection growing radius

# ICFIT parameters used for fitting
define	IC_NFIT		Memi[$1+10]	# Number of fit points
define	IC_NREJECT	Memi[$1+11]	# Number of rejected points
define	IC_RG		Memi[$1+12]	# Pointer for ranges
define	IC_XFIT		Memi[$1+13]	# Pointer to ordinates of fit points
define	IC_YFIT		Memi[$1+14]	# Pointer to abscissas of fit points
define	IC_WTSFIT	Memi[$1+15]	# Pointer to weights of fit points
define	IC_REJPTS	Memi[$1+16]	# Pointer to rejected points
define	IC_METHOD	Memi[$1+17]	# Minimization method

# ICFIT parameters used for interactive graphics
define	IC_NEWX		Memi[$1+18]	# New x fit points?
define	IC_NEWY		Memi[$1+19]	# New y points?
define	IC_NEWWTS	Memi[$1+20]	# New weights?
define	IC_NEWFUNCTION	Memi[$1+21]	# New fitting function?
define	IC_OVERPLOT	Memi[$1+22]	# Overplot next plot?
define	IC_FITERROR	Memi[$1+23]	# Error in fit
define	IC_LABELS	Memi[$1+24+$2-1]# Graph axis labels
define	IC_UNITS	Memi[$1+26+$2-1]# Graph axis units
define	IC_HELP		Memi[$1+28]	# Pointer to help file name
define	IC_COMP		Memi[$1+29]	# Plot individual components ?
define	IC_XAXIS	Memi[$1+30]	# Type of x axis
define	IC_YAXIS	Memi[$1+31]	# Type of y axis
define	IC_MAG0		Memr[$1+32]	# Constant for magnitude y scale
define	IC_MINY		Memr[$1+33]	# Minimum positive y value.
define	IC_LTYPE	Memi[$1+34]	# Continuous line or markers.

# ICFIT key definitions
define	IC_GKEY		Memi[$1+35]			# Graph key
define	IC_AXES		Memi[$1+36+($2-1)*2+$3-1]	# Graph axis codes

# Default help file and prompt
define	IC_DEFHELP	"stsdas$pkg/analysis/fitting/ncfit/icgfit.key"
define	IC_PROMPT	"non-linear icfit cursor options"

# Non-overlapping names with original icfit.
include		"names.h"

