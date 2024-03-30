# GRATING.H -- Variables used to describe grating dispersion

define	LEN_G		14
define	GRATING_STR	"m1,m2,size,f,gy,beta_y,sigma_y,delta_y,theta_y,gx,beta_x,\
sigma_x,delta_x,theta_x"

define	M1	1		# first spectral order
define	M2	2		# last spectral order
define	SIZE	3		# detector size
define	F	4		# focal length
define	GY	5		# grating constant
define	BETA_Y	6		# echelle blaze angle
define	SIGMA_Y	7		# half angle btw collimator and disperser
				# parallel to ruling
define	DELTA_Y	8		# half angle btw collimator and disperser
				# perpendicular to ruling
define	THETA_Y	9		# grating scan angle
define	GX	10		# cross disperser grating constant
define	BETA_X	11		# cross disperser blaze angle
define	SIGMA_X	12		# half angle btw collimator and disperser
				# parallel to ruling
define	DELTA_X	13		# half angle btw collimator and disperser
				# perpendicular to ruling
define	THETA_X	14		# grating scan angle

# Macros for common code dealing wiyh both dispersion directions

define	Y_DIR		1
define	X_DIR		2

define	GCON		$1[5*$2]
define	BETA		$1[5*$2+1]
define	SIGMA		$1[5*$2+2]
define	DELTA		$1[5*$2+3]
define	THETA		$1[5*$2+4]
