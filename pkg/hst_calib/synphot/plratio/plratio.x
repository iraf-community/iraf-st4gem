include "../limit.h"

#* HISTORY *
#* B.Simon	14-Jul-94	original

# PLRATIO -- Plot the ratio between spectra and spectrophotometry files

procedure plratio ()

#--
pointer	obsmode		# observation mode
pointer	spectrum	# target spectrum
pointer	form		# output form
pointer	spfile		# spectrophotometry data
pointer	pfile		# photometry data
pointer	vzero		# variable list
real	inlimit[4]	# user specified plot limits
bool	append		# append to existing plot?
pointer	ltype		# line type
pointer	device		# graphics device
pointer	wavetab		# wavelength table
pointer	grtbl		# graph table
pointer	cmptbl		# component lookup table
real	area		# telescope area

bool	ylog
int	gmode, nstat, pstat, ntgt
pointer	sp, gp, dspec, dsphot, dphot, dratio, tgtlist
real	chisq, bias, rms, pchisq, pbias, prms, outlimit[4]

data	ylog       / false /
string	title      "SYNPHOT.PLRATIO"
string	obslabel   "OBSMODE ="
string	speclabel  "REFSPEC ="
string	vzlabel    "$0 ="
string	splabel    "SPHOT ="
string	photlabel  "PHOT ="
string	rmslabel   "RMS ="
string	biaslabel  "BIAS ="
string	chilabel   "CHISQ ="
string  numlabel   "NUM PTS ="
string	prmslabel  "PHOT RMS ="
string	pbiaslabel "PHOT BIAS ="
string	pchilabel  "PHOT CHISQ ="

bool	clgetb()
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
	call salloc (ltype, SZ_FNAME, TY_CHAR)
	call salloc (device, SZ_FNAME, TY_CHAR)
	call salloc (wavetab, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("obsmode", Memc[obsmode], SZ_FNAME)
	call clgnone ("spectrum", Memc[spectrum], SZ_FNAME)
	call clgstr ("form", Memc[form], SZ_FNAME)
	call clgnone ("spfile", Memc[spfile], SZ_FNAME)
	call clgnone ("pfile", Memc[pfile], SZ_FNAME)
	call clgstr ("vzero", Memc[vzero], SZ_FNAME)

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

	# Compute ratio of spectrophotometry files to spectra

	call getdspec (Memc[obsmode], Memc[spectrum], Memc[form], 
		       Memc [vzero], Memc[wavetab], Memc[grtbl], 
		       Memc[cmptbl], dspec)

	call getdsphot (Memc[spfile], Memc[form], Memc[wavetab], dsphot)

	call getdratio (Memc[wavetab], Memc[form], dspec, dsphot, dratio)

	call statratio (Memc[form], dspec, dsphot, dratio, 
			chisq, bias, rms, nstat)

	# Compute ratio of photometry data to spectra

	call getdphot (Memc[pfile], Memc[form], Memc[wavetab], Memc[grtbl], 
		       Memc[cmptbl], dphot)

	call tgtmatch (Memc[obsmode], Memc[spectrum], Memc [vzero], 
		       Memc[pfile], tgtlist, ntgt)

	call statphot (Memc[form], tgtlist, ntgt, dspec, dphot, 
		       pchisq, pbias, prms, pstat)

	call photratio (Memc[form], tgtlist, ntgt, dspec, dphot)

	# Open graphics device, compute plot limits, and plot frame

	if (append) {
	    gp = gopen (Memc[device], APPEND, STDGRAPH)

	} else {
	    gp = gopen (Memc[device], NEW_FILE, STDGRAPH)

	    call limdratio (dratio, inlimit, outlimit)
	    call limdphot (dphot, Memc[form], inlimit, outlimit)

	    call boxratio (gp, inlimit, outlimit, Memc[form], title)
	}

	# Plot data

	call setmarker (gp, Memc[ltype])
	call pltdratio (gp, Memc[form], dratio)
	call pltdphot (gp, dphot)

	# Label plot

	call labelt (gp, obslabel, Memc[obsmode])
	call labelt (gp, speclabel, Memc[spectrum])
	call labelt (gp, vzlabel, Memc[vzero])
	call labelt (gp, splabel, Memc[spfile])
	call labelt (gp, photlabel, Memc[pfile])

	call labelr (gp, rmslabel, rms)
	call labelr (gp, biaslabel, bias)
	call labelr (gp, chilabel, chisq)
	call labeli (gp, numlabel, nstat)

	call labelr (gp, prmslabel, prms)
	call labelr (gp, pbiaslabel, pbias)
	call labelr (gp, pchilabel, pchisq)
	call labeli (gp, numlabel, pstat)

	# Release memory

	call clsdspec (dspec)
	call clsdsphot (dsphot)
	call clsdphot (dphot)
	call clsdratio (dratio)

	if (tgtlist != NULL)
	    call mfree (tgtlist, TY_INT)

	call gclose (gp)
	call sfree (sp)

end
