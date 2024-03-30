#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	13-Apr-95	added size check
#* B.Simon	03-Jun-99	added keyword existence check

# PIXSIZE -- Compute pixel size of image

procedure pixsize (im, size)

pointer	im		# i: image descriptor
real	size		# o: pixel size
#--
real	r[2], w[2], cd[2,2]
pointer	mw, sp, type1, type2

string	badaxis  "|pixels|"

int	imaccf(), strdic()
pointer	mw_openim()

begin
	# Allocate memory for temporary strings
	call smark (sp)
	call salloc (type1, SZ_FNAME, TY_CHAR)
	call salloc (type2, SZ_FNAME, TY_CHAR)

	# Calculate pixel size from cdmatrix

	mw = mw_openim (im)
	call mw_gwtermr (mw, r, w, cd, 2)

	size = abs (cd[1,1] * cd[2,2]  - cd[1,2] * cd[2,1])
	size = sqrt (size)

	call mw_close (mw)

	# Check for bad axis types

	if (imaccf (im, "CTYPE1") == NO || imaccf (im, "CTYPE2") == NO) {
	    size = 0.0

	} else {
	    call imgstr (im, "CTYPE1", Memc[type1], SZ_FNAME)
	    call imgstr (im, "CTYPE2", Memc[type2], SZ_FNAME)
	    call strfix (Memc[type1])
	    call strfix (Memc[type2])

	    if (strdic (Memc[type1], Memc[type1], SZ_FNAME, badaxis) != 0)
		size = 0.0

	    if (strdic (Memc[type2], Memc[type2], SZ_FNAME, badaxis) != 0)
		size = 0.0
	}

	call sfree (sp)
end

