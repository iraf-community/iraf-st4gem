# LIBSYNPHOT.H -- Internal constants used by the synphot library

define	HORNE		NO		# old synphot compatibility mode

define	NO_COLUMN	0		# no column number

define	MISC_ERR	1		# miscellaneous error
define	SYN_ERR		2		# expression syntax error
define	COL_ERR		3		# column not found error

define  GRF_COMPID      1               # index to optical component col in graph tab
define  GRF_KEYWRD      2               # index to keyword col in graph tab
define  GRF_INNODE      3               # index to innode col in graph tab
define  GRF_OUTNODE     4               # index to outnode col in graph tab
define  GRF_THMLID      5               # index to thermal component col in graph tb

define	MAXCACHE	12		# Maximum number of cached tables
define	MAXLIST		100		# Maximum number of components in list
define	MAXPARAM	3		# Max number of parameters / component

define	SZ_KEYWRD	32		# Length of keyword
define	SZ_COMPID      	18		# Length of component name

define	PCH	        '#'		# Leading character of parameter

define	LENWAVE		10000		# Length of default wavelength set

define	PI		3.14159265	# Mysterious math constant
define	H		6.62620e-27	# Planck's constant
define	C		2.997925e18	# Speed of light
define	KBOLTZ		1.3806e-16	# Boltzman's constant
define	RSUN 		6.9599e10	# Radius of sun
define	PC 		3.085678e18	# Parsec

define	RADIAN		(RSUN/PC/1000.)
define	RENORM		(PI*RADIAN*RADIAN)

define  ABZERO		(-48.60)  	# zero point of the AB magnitude system
define  STZERO		(-21.10)  	# zero point of the ST magnitude system

define	FORMSTR		"photlam,counts,flam,fnu,photnu,jy,mjy,\
abmag,stmag,vegamag,obmag"		# list of valid flux units

define	VEGA		"st4gem$lib/synphot/vega.dat"

