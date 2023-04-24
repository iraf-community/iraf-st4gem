define	ARGEMAX	      88.D+0	# Maximum argument of exp()

# NL_2DGAUSS -- Compute 2-d Gaussian.
#
# The 2.7725 factor relates to the Gaussian being expressed as a
# function of FWHM instead of standard deviation:
#
#  2 * stddev = FWHM / sqrt (2 * ln (2))

procedure nl_2dgauss (x, y, z, ndata, amp, x0, y0, fwhm, ellip, teta)

real	x[ARB]		# i: Independent x variable array
real	y[ARB]		# i: Independent y variable array
real	z[ARB]		# io: Output array with Gaussian
int	ndata		# i: Number of data points.
real	amp		# i: Intensity
real	x0		# i: X coordinate of center
real	y0		# i: Y coordinate of center
real	fwhm		# i: FWHM
real	ellip		# i: ellipticity
real	teta		# i: position angle (radians)

real	x1, y1
double	arg, radius
int	i
bool	flag

begin
	if (ellip == 0.)
	    flag = false
	else
	    flag = true
	do i = 1, ndata {
	    if (flag) {
	        x1 = (x[i] - x0) * cos(teta) + (y[i] - y0) * sin(teta)
	        y1 = (y[i] - y0) * cos(teta) - (x[i] - x0) * sin(teta)
	        radius = double (sqrt (x1**2 * (1. - ellip)**2 + 
				       y1**2) / (1. - ellip))
	    } else {
	        radius = double (sqrt ((x[i] - x0)**2 + (y[i] - y0)**2))
	    }
	    arg = -2.7725d0 * (radius / fwhm) ** 2
	    arg = min(arg, ARGEMAX)
	    z[i] = amp * exp (arg)
	}
end
                                 
     
