include	<math.h>

#* HISTORY *
#* B.Simon	17-Feb-95	original

# RADCALC -- Calculate the radius for an axially  symmetric function

procedure radcalc (x, y, r)

real	x	# i: x position
real	y	# i: y position
real	ar	# i: axial ratio
real	pa	# i: position angle
real	r	# o: radius
#--
real	a, b, c
real	rpa, cpa, spa, t

string	badax  "Illegal value of axial ratio"

begin

	r = sqrt (a * x ** 2 + b * x * y + c * y ** 2)
	return

	entry initrad (ar, pa)

	if (ar <= 0.0 )
	    call printerr_real (badax, ar)

	rpa = DEGTORAD (pa)

 	cpa = cos (rpa)
	spa = sin (rpa)
	t = 1.0 / (ar * ar)

	a = cpa * cpa + spa * spa * t
	b = 2.0 * (1.0 - t) * cpa * spa
	c = spa * spa + cpa * cpa * t

	return
end
