# This include file contains definitions used in tasks that handle 
# linear curve fitting.


# Name of straight power-series polynomial mode.

define	PSPOLY		"poly"


# The following definition MUST be the same as in 
# iraf$pkg/xtools/icfit/icparams.x !!!!

define	FUNCTIONS	"|chebyshev|legendre|spline3|spline1|aaaa"

define	CV_CHEB		1
define	CV_LEG		2
define	CV_SP3		3
define	CV_SP1		4
