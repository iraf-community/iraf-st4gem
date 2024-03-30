include "../limit.h"

#* HISTORY *
#* B.Simon	16-Jun-94	original

# PLSPEC -- Plot spectra

procedure plspec ()

#--
pointer	obsmode		# observation mode
pointer	spectrum	# target spectrum
pointer	form		# output form
pointer	vzero		# variable list
pointer	spfile		# spectrophotometry data
pointer	pfile		# photometry data
pointer	errtyp		# error type
real	inlimit[4]	# user specified plot limits
bool	append		# append to existing plot?
pointer	ltype		# line type
pointer	device		# graphics device
pointer	wavetab		# wavelength table
pointer	grtbl		# graph table
pointer	cmptbl		# component lookup table
real	area		# telescope area

bool	ylog
int	gmode
pointer	sp, gp, dspec, dsphot, dphot
real	outlimit[4]

data	ylog       / false /
string	title      "SYNPHOT.PLSPEC"
string	obslabel   "OBSMODE ="
string	speclabel  "REFSPEC ="
string	vzlabel    "$0 ="
string	splabel    "SPHOT ="
string	photlabel  "PHOT ="
string	alabel     "AREA (cm^2) ="

bool	clgetb()
int	is_count()
pointer	gopen()
real	clgetr()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (obsmode, SZ_FNAME, TY_CHAR)
	call salloc (spectrum, SZ_FNAME, TY_CHAR)
	call salloc (form, SZ_FNAME, TY_CHAR)
	call salloc (vzero, SZ_FNAME, TY_CHAR)
	call salloc (spfile, SZ_FNAME, TY_CHAR)
	call salloc (pfile, SZ_FNAME, TY_CHAR)
	call salloc (errtyp, SZ_FNAME, TY_CHAR)
	call salloc (ltype, SZ_FNAME, TY_CHAR)
	call salloc (device, SZ_FNAME, TY_CHAR)
	call salloc (wavetab, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("obsmode", Memc[obsmode], SZ_FNAME)
	call clgnone ("spectrum", Memc[spectrum], SZ_FNAME)
	call clgstr ("form", Memc[form], SZ_FNAME)
	call clgstr ("vzero", Memc[vzero], SZ_FNAME)
	call clgnone ("spfile", Memc[spfile], SZ_FNAME)
	call clgnone ("pfile", Memc[pfile], SZ_FNAME)
	call clgstr ("errtyp", Memc[errtyp], SZ_FNAME)

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

	call clgstr ("ltype", Memc[ltype], SZ_FNAME)
	call clgstr ("device", Memc[device], SZ_FNAME)
	call clgnone ("wavetab", Memc[wavetab], SZ_FNAME)
	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)

	area = clgetr ("area")
	call put_hstarea (area)

	# Compute quantities to be plotted

	call getdspec (Memc[obsmode], Memc[spectrum], Memc[form], 
		       Memc [vzero], Memc[wavetab], Memc[grtbl], 
		       Memc[cmptbl], dspec)

	call getdsphot (Memc[spfile], Memc[form], 
			Memc[wavetab], dsphot)

	call getdphot (Memc[pfile], Memc[form], Memc[wavetab], 
		       Memc[grtbl], Memc[cmptbl], dphot)

	# Open graphics device, compute plot limits, and plot frame

	if (append) {
	    gp = gopen (Memc[device], APPEND, STDGRAPH)

	} else {
	    gp = gopen (Memc[device], NEW_FILE, STDGRAPH)

	    call limdspec (dspec, Memc[form], inlimit, outlimit)
	    call limdsphot (dsphot, inlimit, outlimit)
	    call limdphot (dphot, Memc[form], inlimit, outlimit)

	    call pltbox (gp, ylog, inlimit, outlimit, Memc[form], title)
	}

	# Plot data

	call setmarker (gp, Memc[ltype])

	call pltdspec (gp, Memc[form], dspec)
	call pltdsphot (gp, dsphot, Memc[errtyp])
	call pltdphot (gp, dphot)

	# Label plot

	call labelt (gp, obslabel, Memc[obsmode])
	call labelt (gp, speclabel, Memc[spectrum])
	call labelt (gp, vzlabel, Memc[vzero])
	call labelt (gp, splabel, Memc[spfile])
	call labelt (gp, photlabel, Memc[pfile])

	if (is_count (Memc[form]) == YES)
	    call labelr (gp, alabel, area)

	# Release memory

	call clsdspec (dspec)
	call clsdsphot (dsphot)
	call clsdphot (dphot)

	call gclose (gp)
	call sfree (sp)

end
