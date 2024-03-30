# USRLIM -- Set plot limits

procedure usrlim( xmin, xmax, ymin, ymax )

real	xmin	# o: minimum x value
real	xmax	# o: maximum x value
real	ymin	# o: minimum y value
real	ymax	# o: maximum y value
#--
real	clgetr()

begin

	xmin = clgetr( "left" )
	xmax = clgetr( "right" )

	ymin = clgetr( "bottom" )
	ymax = clgetr( "top" )

end
