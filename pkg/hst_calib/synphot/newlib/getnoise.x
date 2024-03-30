include	<math.h>

#* HISTORY *
#* B.Simon	30-Mar-95	original

# GETNOISE -- Compute poisson random noise from its mean

procedure getnoise (seed, nmean, mean, noise)

long	seed		# i: random number seed
int	nmean		# i: length of mean array
real	mean[ARB]	# i: mean noise
real	noise[ARB]	# o: random noise
#--
int	im
real	omean, val, g, t, sqmn, logmn, y

real	urand(), log_gamma()

begin
	# The algorithm is based on POIDEV in section 7.3 of
	# Numerical Recipies

	omean = -1.0
	do im = 1, nmean {
	    if (mean[im] <= 0.0) {
		# Set deviate to zero for zero or negative mean
		val = 0.0

	    } else if (mean[im] <= 12.0) {
		# Compute deviate by direct method 
		# for small values of the mean

		if (mean[im] != omean) {
		    omean = mean[im]
		    g = exp (- omean)
		}

		t = 1.0
		val = -1.0
		repeat {
		    val = val + 1.0
		    t =  t * urand (seed)
		} until (t <= g)

	    } else {
		# Compute deviate by rejection method
		# for large values of the mean

		# Compute constants used in method

		if (mean[im] != omean) {
		    omean = mean[im]
		    sqmn = sqrt (2.0 * omean)
		    logmn = alog (omean)
		    g = omean * logmn - log_gamma (omean + 1.0)
		}

		repeat {
		    # Lorentzian comparison function

		    repeat {
			y = tan (PI * urand (seed))
			val = y * sqmn + omean
		    } until (val >= 0.0)

		    # Compare to poisson distribution

		    val = int (val)
		    t = 0.9 * (1.0 + y ** 2) *
			exp (val * logmn - log_gamma (val + 1.0) - g)

		} until (urand (seed) <= t)
	    }

	    noise[im] = val
	}

end

# LOG_GAMMA -- Log of gamma function by Lanczos approximation

real procedure log_gamma (val)

real	val		# i: function argument
#--
double	half, one, gfac, rtwopi, coef[6]
double	arg, factor, term
int	i

data	half / 0.5d0 /
data	one  / 1.0d0 /
data 	gfac / 5.5d0 /

data	rtwopi / 2.50662827465d0 /

data	coef / 76.18009173d0, -86.50532033d0, 24.01409822d0, 
	       -1.231739516d0, .120858003d-2, -.536382d-5 /

begin
        # Lanczos approximation for gamma function is from GAMMLN
        # in section 6.1 of Numerical Recipies

	# Compute leading factor 

	arg = val - one
	factor = arg + gfac
	factor = (arg + half) * log (factor) - factor

	# Compute series expansion

	term = one
	do i = 1, 6 {
	    arg = arg + one
	    term = term + coef[i] / arg
	}

	return (factor + log (rtwopi * term))
end
