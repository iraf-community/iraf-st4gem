include	<mach.h>

# IS_LINEAR -- Check to see whether wavelength array is spaced linearly

bool procedure is_linear (wave, nwave)

real	wave[ARB]	# i: wavelength array
int	nwave		# i: number of wavelengths
#--
bool	linear
int	iwave
real	tol, dif2

begin

	linear = true
	tol = 10.0 * EPSILONR

	do iwave = 1, nwave-2 {
	    dif2 = wave[iwave+2] - 2.0 * wave[iwave+1] + wave[iwave]
	    if (abs (dif2 / wave[iwave+1]) > tol) {
		linear = false
		break
	    }
	}

	return (linear)
end
