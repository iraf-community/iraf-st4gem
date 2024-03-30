include <math.h>

#* HISTORY *
#* B.Simon	14-Apr-95	original

# ECLIPTIC -- Compute the ecliptic coordinates of a target

procedure ecliptic (ra, dec, lambda, beta)

double	ra		# i: target right ascension (degrees)
double	dec		# i: target declination (degrees)
double	lambda		# o: eccliptic longitude
double	beta		# o: eccliptic lattiutde
#--
double	eps, alpha, delta

data	eps / 23.4392911 /	# obliquity of ecliptic for J2000

begin
	# Formulas are from "Astronomical Algorithms" by Jean Meeus
	# chapter 12

	# Convert angles from degrees to radians

	alpha = DEGTORAD(ra)
	delta = DEGTORAD(dec)
	eps = DEGTORAD (eps)

	# Compute ecliptic longitude and latitude

	lambda = atan2 (sin (alpha) * cos (eps) + tan (delta) * sin (eps),
			cos (alpha))

	beta = asin (sin (delta) * cos (eps) - 
		     cos (delta) * sin (eps) * sin (alpha))

	# Convert from radians back to degrees

	lambda = RADTODEG (lambda)
	beta = RADTODEG (beta)

end
