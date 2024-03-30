#* HISTORY *
#* B.Simon	17-Feb-95	original

# PROFINTERP -- Interpolate a profile function on an array

procedure profinterp (func, nsub, nx, ny, array)

real	func		# i: function to compute profile
int	nsub		# i: number of pixel subdivisions
int	nx		# i: array size
int	ny		# i: array size
real	array[nx,ny]	# o: output array
#--
int	ix, iy
pointer	sp, xgrid, ygrid
real	scale, total, xlo, xhi, ylo, yhi, sum

extern	func

begin
	call smark (sp)
	call salloc (xgrid, nx+1, TY_REAL)
	call salloc (ygrid, ny+1, TY_REAL)

	scale =  1.0 / nsub
	call setgrid (scale, Memr[xgrid], nx+1)
	call setgrid (scale, Memr[ygrid], ny+1)

	total = 0.0
	do iy = 1, ny {
	    ylo = Memr[ygrid+iy-1]
	    yhi = Memr[ygrid+iy]
	    sum = INDEFR

	    do ix = 1, nx {
		xlo = Memr[xgrid+ix-1]
		xhi = Memr[xgrid+ix]

		call simpson2 (xlo, xhi, ylo, yhi, func, sum, array[ix,iy])
		total = total + array[ix,iy]
	    }
	}

	if (total > 0.0) {
	    total = 1.0 / total
	    call amulkr (array, total, array, nx*ny)

	} else {
	    ix = (nx + 1) / 2
	    iy = (ny + 1) / 2
	    array[ix,iy] = 1.0
	}

	call sfree (sp)
end
