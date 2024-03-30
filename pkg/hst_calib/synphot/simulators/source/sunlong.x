include	<math.h>

#* HISTORY *
#* B.Simon	14-Apr-95	original

# SUNLONG -- Compute the sun's ecliptic longitude in degrees

procedure sunlong (jd, lambda)

double	jd		# i: julian date
double	lambda		# o: sun's ecliptic longitude
#--
double	t, m, c, lz, omega

begin
	# Formulas are from "Astronomical Algorithms" by Jean Meeus
	# chapter 24

	# Compute centuries from 2000

	t = (jd - 2451545.0) / 36525.0

	# Compute mean anomaly

	m = 357.52910 + 35999.05030 * t - 0.0001559 * t * t -
	    0.00000048 * t * t * t
	m = mod (m, 360.0d0)
	m = DEGTORAD(m)

	# Compute equation of center

	c = (1.914600 - 0.004817 * t - 0.000014 * t * t) * sin (m) +
	    (0.019993 - 0.000101 * t) * sin (2. * m) +
	     0.000290 * sin (3. * m)

	# Compute mean longitude, add to equation of center to give
	# true longitude

	lz = 280.46645 + 36000.76983 * t + 0.0003032 * t * t
	lz = mod (lz, 360.0d0)
	lambda = lz + c

	# Compute apparent longitude

	omega = 125.04 - 1934.136 * t
	lambda = lambda - 0.00569 - 0.00478 * sin (DEGTORAD(omega))

	# Convert from equinox of date to J2000

	lambda = lambda - 1.397 * t

end
