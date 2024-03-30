include "../limit.h"

#* HISTORY *
#* B.Simon	21-Jul-94	original

# PLTRANS -- Plot a transformation diagram

procedure pltrans ()

#--
pointer	spectrum	# target spectrum
pointer	xmode		# x axis mode
pointer	ymode		# y axis mode
pointer	xform		# x axis form
pointer	yform		# y axis form
pointer	vzero		# variable list
pointer	input		# input table containing photometric data
pointer	output		# output table with photometric data
real	inlimit[4]	# user specified plot limits
bool	append		# append to existing plot?
pointer	mktype		# marker type
pointer	device		# graphics device
pointer	wavetab		# wavelength table
pointer	grtbl		# graph table
pointer	cmptbl		# component lookup table
real	area		# telescope area

int	gmode
pointer	sp, gp, dtrans
real	outlimit[4]

string	title      "SYNPHOT.PLTRANS"

bool	clgetb()
pointer	gopen()
real	clgetr()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (spectrum, SZ_FNAME, TY_CHAR)
	call salloc (xmode, SZ_FNAME, TY_CHAR)
	call salloc (ymode, SZ_FNAME, TY_CHAR)
	call salloc (xform, SZ_FNAME, TY_CHAR)
	call salloc (yform, SZ_FNAME, TY_CHAR)
	call salloc (vzero, SZ_FNAME, TY_CHAR)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (mktype, SZ_FNAME, TY_CHAR)
	call salloc (device, SZ_FNAME, TY_CHAR)
	call salloc (wavetab, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("spectrum", Memc[spectrum], SZ_FNAME)
	call clgnone ("xmode", Memc[xmode], SZ_FNAME)
	call clgnone ("ymode", Memc[ymode], SZ_FNAME)
	call clgstr ("xform", Memc[xform], SZ_FNAME)
	call clgstr ("yform", Memc[yform], SZ_FNAME)
	call clgstr ("vzero", Memc[vzero], SZ_FNAME)
	call clgnone ("input", Memc[input], SZ_FNAME)
	call clgnone ("output", Memc[output], SZ_FNAME)

	inlimit[LEFT] = clgetr ("left")
	inlimit[RIGHT] = clgetr ("right")
	inlimit[BOTTOM] = clgetr ("bottom")
	inlimit[TOP] = clgetr ("top")
	call amovkr (INDEFR, outlimit, 4)

	append = clgetb ("append")

	if (append) {
	    gmode = APPEND
	} else {
	    gmode = NEW_FILE
	}

	call clgstr ("mktype", Memc[mktype], SZ_FNAME)
	call clgstr ("device", Memc[device], SZ_FNAME)
	call clgnone ("wavetab", Memc[wavetab], SZ_FNAME)
	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)

	area = clgetr ("area")
	call put_hstarea (area)

	call chkform (Memc[xform])
	call chkform (Memc[yform])

	# Compute transformation diagram data

	call inpdtrans (Memc[input], Memc[xmode], Memc[ymode], Memc[xform], 
		        Memc[yform], SZ_FNAME, dtrans)

	call getdtrans (Memc[spectrum], Memc[xmode], Memc[ymode], Memc[xform], 
		        Memc[yform], Memc [vzero], Memc[wavetab], Memc[grtbl], 
			Memc[cmptbl], dtrans)

	call wrtdtrans (Memc[output], Memc[xmode], Memc[ymode], Memc[xform], 
		        Memc[yform], append, dtrans)

	# Open graphics device, compute plot limits, and plot frame

	if (append) {
	    gp = gopen (Memc[device], APPEND, STDGRAPH)

	} else {
	    gp = gopen (Memc[device], NEW_FILE, STDGRAPH)

	    call limdtrans (dtrans, inlimit, outlimit)

	    call boxtrans (gp, inlimit, outlimit, Memc[xmode], Memc[ymode], 
			   Memc[xform], Memc[yform], title)
	}

	# Plot data

	call pltdtrans (gp, Memc[mktype], dtrans)

	# Release memory

	call clsdtrans (dtrans)

	call gclose (gp)
	call sfree (sp)

end
