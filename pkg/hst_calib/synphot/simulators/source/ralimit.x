define	CYCLE		360.0

#* HISTORY *
#* B.Simon	17-Feb-95	original

# RALIMIT -- Put object RA within 180 degree limits of target RA

procedure ralimit (target, object)

double	target		# i: target RA in degrees
double	object		# u: object RA in degrees
#--
double	hi, lo

begin
	hi = target + CYCLE / 2.0
	lo = target - CYCLE / 2.0

	while (object > hi)
	    object = object - CYCLE

	while (object < lo)
	    object = object + CYCLE
end
