#* HISTORY *
#* B.Simon	10-Nov-95	original

# POLYAREA -- Area of a polygon

real procedure polyarea (npoly, poly)

int	npoly		# i: number of polygon vertices
real	poly[2,npoly]	# i: polygon vertices
#--
int	ipoly
real	area, x1, x2, y1, y2

begin
	area = 0.0

	# First vertex is repeated as last vertex, so skip it

	do ipoly = 2, npoly - 2 {
	    # Divide polygon into triangles using first vertex as origin

	    x1 = poly[1,ipoly] - poly[1,1]
	    x2 = poly[1,ipoly+1] - poly[1,1]

	    y1 = poly[2,ipoly] - poly[2,1]
	    y2 = poly[2,ipoly+1] - poly[2,1]

	    # Compute area of each triangle using cross product
	    # Polygon area is sum of triangular areas

	    area = area + 0.5 * abs (x1 * y2 - x2 * y1)
	}

	return (area)
end
