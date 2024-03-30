include	<gset.h>
include	<mach.h>
include	"../limit.h"
define	FORMSTR   "photlam,counts,flam,fnu,photnu,jy,mjy,abmag,stmag,\
vegamag,obmag"

#* HISTORY *
#* B.Simon	09-Jun-94	Original

# PLTBOX -- Plot and label the axes and plot the title

procedure pltbox (gp, ylog, inlimit, outlimit, form, title)

pointer	gp		# i: graphics descriptor
bool	ylog		# i: Semi-log plot?
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# i: calculated plot limits
char	form[ARB]	# i: spectral form
char	title[ARB]	# i: plot title
#--
char 	formlab[31,11]
int	isband, ftype
pointer	sp, ylabel
real	dif, temp, limit[4]

int	word_match(), is_magunit()

string	formlist  FORMSTR
string	xlabel    "WAVELENGTH (A)"
string	noextent  "No vertical extent in plot"

string	formlab1  "NLAM (PHOTONS cm-2 s-1 A-1)"
string	formlab2  "DETECTED COUNTS (Area-1 s-1)"
string	formlab3  "FLAM (ergs cm-2 s-1 A-1)"
string	formlab4  "FNU (ergs cm-2 s-1 Hz-1)"
string	formlab5  "NNU (PHOTONS cm-2 s-1 Hz-1)"
string	formlab6  "FNU (Jy)"
string	formlab7  "FNU (mJy)"
string	formlab8  "ABNU (MAGS)"
string	formlab9  "STLAM (MAGS)"
string	formlabz  "VEGAMAG"
string	formlaby  "OBMAG"

equivalence (formlab[1,1], formlab1)
equivalence (formlab[1,2], formlab2)
equivalence (formlab[1,3], formlab3)
equivalence (formlab[1,4], formlab4)
equivalence (formlab[1,5], formlab5)
equivalence (formlab[1,6], formlab6)
equivalence (formlab[1,7], formlab7)
equivalence (formlab[1,8], formlab8)
equivalence (formlab[1,9], formlab9)
equivalence (formlab[1,10], formlabz)
equivalence (formlab[1,11], formlaby)

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (ylabel, SZ_FNAME, TY_CHAR)

	# Get ylabel from spectral form

	call strcpy (form, Memc[ylabel], SZ_FNAME)
	call strfix (Memc[ylabel])

	if (Memc[ylabel] == EOS) {
	    isband = YES
	    call strcpy ("PASSBAND", Memc[ylabel], SZ_FNAME)

	} else {
	    isband = NO
	    ftype = word_match (Memc[ylabel], formlist)
	    if (ftype == 0) {
		call strupr (Memc[ylabel])
	    } else {
		call strcpy (formlab[1,ftype], Memc[ylabel], SZ_FNAME)
	    }
	}

	# Adjust plot limits

	call amovr (outlimit, limit, 4)

	dif = limit[TOP] - limit[BOTTOM]
	if (dif == 0.0)
	    call printerr_real (noextent, limit[BOTTOM])

	if (IS_INDEFR(inlimit[BOTTOM]))
	    limit[BOTTOM] = limit[BOTTOM] - 0.05 * dif

	if (IS_INDEFR(inlimit[TOP]))
	    limit[TOP] = limit[TOP] + 0.05 * dif

	# Swap limits for magnitude units

	if (is_magunit (form) == YES) {
	    temp = limit[TOP]
	    limit[TOP] = limit[BOTTOM]
	    limit[BOTTOM] = temp

	} else if (ylog) {
	    limit[BOTTOM] = max (10.0 / MAX_REAL, limit[BOTTOM])
	}

	# Set limits and plot axes

	call gswind (gp, limit[LEFT], limit[RIGHT], limit[BOTTOM], limit[TOP])

	call gseti (gp, G_TXQUALITY, GT_HIGH)
	if (ylog)
	    call gseti (gp, G_YTRAN, GW_LOG)

	call glabax (gp, title, xlabel, Memc[ylabel])

	# Plot zero line for spectra if form is not magnitude

	if (isband == NO && is_magunit (form) == NO) {
	    call gamove (gp, limit[LEFT], 0.0)
	    call gadraw (gp, limit[RIGHT], 0.0)
	}

	# Unset cursor position for label

	call labelu
	call sfree (sp)
end
