#* HISTORY
#* B.Simon	22-Apr-97	original

# STEPOFFSET -- Calculate offset needed to step coefficients 
#
# The step coeffiecients need to align with the start of the wavelength array.
# The two may be different if the wavelength array is truncated beacuse of
# the value selected for the central wavelength.

int procedure stepoffset (startwave, a, b, c)

real	startwave	# i: wavelength that coefficients must align with
real	a		# i: quadratic term in step coefficients
real	b		# i: linear term in step coefficients
real 	c		# i: constant term in step coefficients
#--
int	offset
real	diff, x1, x2

begin
	# Shift the constant term to 
	diff = startwave - c

	if (a == 0.0) {
	    # Solution for linear equation. The 0.5 rounds the result
	    # to the nearest integer

	    offset = diff / b + 0.5 

	} else {
	    # Alternate form of quadratic equation that avoids ill conditioning

	    x1 = - (b + sqrt (b * b + 4.0 * a * diff)) / (2.0 * a)
	    x2 = - diff / (a * x1)

	    # The solution closest to zero is the proper solution
	    # Adding the 0.5 rounds to the nearest integer

	    if (abs(x1) < abs (x2)) {
		offset = x1 + 0.5
	    } else {
		offset = x2 + 0.5
	    }
	}

	return (offset)
end
