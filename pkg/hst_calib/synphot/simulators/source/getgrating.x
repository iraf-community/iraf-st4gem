#* HISTORY *
#* B.Simon	13-Jul-95	original

include "grating.h"

# GETGRATING -- Read grating parameters from file

procedure getgrating (gratecat, obsmode, grating)

char	gratecat[ARB]	# i: grating catalog
char	obsmode[ARB]	# i: observation mode
real	grating[LEN_G]	# o: grating parameters
#--
int	jrow, icol, ic
pointer	sp, catalog, cname, tp, cp

string	obscol     "OBSMODE"
string	gratecols  GRATING_STR  
string	notfound   "Observation mode not found in grating catalog"

pointer	tbtopn(), word_fetch()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catalog, SZ_FNAME, TY_CHAR)
	call salloc (cname, SZ_FNAME, TY_CHAR)

	# Open aperture dimension  catalog

	call lastfile (gratecat, Memc[catalog], SZ_FNAME)
	tp = tbtopn (Memc[catalog], READ_ONLY, 0)

	# Find table row which matches obsmode

	call findmode (tp, obscol, obsmode, jrow)

	# Copy grating parameters from table

	if (jrow == 0)
	    call printerr_str (notfound, obsmode)

	ic = 1
	icol = 1
	while (word_fetch (gratecols, ic, Memc[cname], SZ_FNAME) > 0) {
	    call syncolptr (tp, Memc[cname], icol+1, cp)
	    call tbegtr (tp, cp, jrow, grating[icol])

	    icol = icol + 1
	}

	# Truncate spectral orders

	grating[M1] = aint (grating[M1])
	grating[M2] = aint (grating[M2])

	if (grating[GY] > 0.0)
	    grating[GY] = 1.0e7 / grating[GY]

	if (grating[GX] > 0.0)
	    grating[GX] = 1.0e7 / grating[GX]

	call tbtclo (tp)
	call sfree (sp)
end

