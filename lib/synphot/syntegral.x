include "libsynphot.h"

# SYNTEGRAL -- Numerical integration of a function

real procedure syntegral (nwave, wave, func)

int	nwave		# i: length of wavelength and function arrays
real	wave[ARB]	# i: wavelength array (independent variable)
real	func[ARB]	# i: tabulated function of wavelength
#--
int	iw
pointer	sp, f2
real	w, wm, wp, sum, dwave

begin
	call smark (sp)
	call salloc (f2, nwave, TY_REAL)

	if (HORNE == NO) {
	    # Piecewise integration of cubic spline. See "Methods of
	    # Numerical Integration", Davis and Rabinowitz, pp. 62-70
	    # Factor of 6 difference in correction term is due to the
	    # different ways of computing the spline coefficients in
	    # this book and "Algorithms" by Sedgewick.

	    call synspline (nwave, wave, func, Memr[f2])

	    sum = 0.0
	    do iw = 1, nwave-1  {
		dwave = wave[iw+1] - wave[iw]
		sum = sum + 0.5 * dwave * (func[iw+1] + func[iw])
		sum = sum - 0.25 * dwave ** 3 * (Memr[f2+iw] + Memr[f2+iw-1])
	    }

	} else {
	    # Keith Horne's midpoint rectangular integration 
	    # provided for backwards compatibility

	    w = wave[1]
	    wp = wave[2]
	    sum = func[1] * (wp - w)

	    do iw = 2, nwave-1 {
		wm = w
		w = wp
		wp = wave[iw+1]
		sum = sum + func[iw] * (wp - wm)
	    }

	    sum = sum + func[nwave] * (wp - w)
	    sum = 0.5 * sum
	}

	call sfree (sp)
	return (sum)
end
