#* HISTORY *
#* B.Simon	17-Feb-95	original

# SETGRID -- Compute an evenly spaced set of points centered on zero

procedure setgrid (scale, grid, ngrid)

real	scale		# i: grid scale
real	grid[ARB]	# o: output grid
int	ngrid		# i: grid length
#--
int	i
real	offset

begin
	offset = -0.5 * (ngrid - 1) * scale
	do i = 1, ngrid
	    grid[i] = (i - 1) * scale + offset

end

