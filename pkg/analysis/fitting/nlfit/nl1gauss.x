define	ARGEMAX	      88.D+0	# Maximum argument of exp()

# NL_1GAUSS -- Adds one Gaussian to z array.
#
# The 2.7725 factor relates to the Gaussian being expressed as a 
# function of FWHM instead of standard deviation:
#
#  2 * stddev = FWHM / sqrt (2 * ln (2))

procedure nl_1gauss (x, z, ndata, amp, center, fwhm)

real	x[ARB]		# i: Independent variable array
real	z[ARB]		# io: Output array with Gaussian added
int	ndata		# i: Number of data points.
real	amp		# i: Intensity
real	center		# i: X value of peak
real	fwhm		# i: FWHM

#--

double	arg
int	i

begin
	if (fwhm <= 0.)
	    return
	else {
	    do i = 1, ndata {
	        arg = -2.7725d0 * ((x[i] - center) / fwhm) ** 2
	        arg = min(arg, ARGEMAX)
	        z[i] = z[i] + amp * exp (arg)
	    }
	}
end
                                 
     
