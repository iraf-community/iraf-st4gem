include "../limit.h"

#* HISTORY *
#* B.Simon	15-Jun-94	original

# PLBAND -- Plot passbands

procedure plband ()

#--
pointer	obsmode		# observation mode
real	inlimit[4]	# user specified plot limits
bool	norm		# normalize passbands to one?
bool	ylog		# take log of y values?
bool	append		# append to existing plot?
pointer	ltype		# line type
pointer	device		# graphics device
pointer	wavetab		# wavelength table
pointer	grtbl		# graph table
pointer	cmptbl		# component lookup table
real	area		# telescope area

int	gmode
pointer	sp, gp, dband
real	outlimit[4]

string	blank     ""
string	title     "SYNPHOT.PLBAND"
string	obslabel  "OBSMODE ="

bool	clgetb()
pointer	gopen()
real	clgetr()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (obsmode, SZ_FNAME, TY_CHAR)
	call salloc (ltype, SZ_FNAME, TY_CHAR)
	call salloc (device, SZ_FNAME, TY_CHAR)
	call salloc (wavetab, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("obsmode", Memc[obsmode], SZ_FNAME)

	inlimit[LEFT] = clgetr ("left")
	inlimit[RIGHT] = clgetr ("right")
	inlimit[BOTTOM] = clgetr ("bottom")
	inlimit[TOP] = clgetr ("top")
	call amovkr (INDEFR, outlimit, 4)

	norm = clgetb ("normalize")
	ylog = clgetb ("ylog")
	append = clgetb ("append")

	if (append) {
	    gmode = APPEND
	} else {
	    gmode = NEW_FILE
	}

	call clgstr ("ltype", Memc[ltype], SZ_FNAME)
	call clgstr ("device", Memc[device], SZ_FNAME)
	call clgnone ("wavetab", Memc[wavetab], SZ_FNAME)
	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)

	area = clgetr ("area")
	call put_hstarea (area)

	# Substitute one for empty obsmode

	if (Memc[obsmode] == EOS)
	    call strcpy ("1.0", Memc[obsmode], SZ_FNAME)

	# Calculate passbands from observation mode

	call getdband (Memc[obsmode], norm, Memc[wavetab], 
		       Memc[grtbl], Memc[cmptbl], dband)

	# Open graphics device, compute plot limits, and plot frame

	if (append) {
	    gp = gopen (Memc[device], APPEND, STDGRAPH)

	} else {
	    gp = gopen (Memc[device], NEW_FILE, STDGRAPH)

	    call limdband (dband, inlimit, outlimit)
	    call pltbox (gp, ylog, inlimit, outlimit, blank, title)
	}

	# Plot passbands

	call setmarker (gp, Memc[ltype])
	call pltdband (gp, ylog, dband)

	call labelt (gp, obslabel, Memc[obsmode])

	# Release memory

	call clsdband (dband)
	call gclose (gp)
	call sfree (sp)

end
