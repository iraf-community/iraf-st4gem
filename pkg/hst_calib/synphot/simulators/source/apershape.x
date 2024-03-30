#* HISTORY *
#* B.Simon	05-Jul-95	original

# APERSHAPE -- Read shape of aperture from catalog

procedure apershape (apercat, obsmode, apshape)

char	apercat[ARB]	# i: catalog of aperture shapes
char	obsmode[ARB]	# i: observation mode
pointer	apshape		# o: aperture shape descriptor
#--
int	jrow
pointer	sp, catalog, func, tp, cp

string	obscol  "OBSMODE"
string	funccol "SHAPE"
string	badmode "No match for obsmode found in aperture shape catalog"

pointer	tbtopn(), breakfunc()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catalog, SZ_FNAME, TY_CHAR)
	call salloc (func, SZ_FNAME, TY_CHAR)

	# Open aperture shape  catalog

	call lastfile (apercat, Memc[catalog], SZ_FNAME)
	tp = tbtopn (Memc[catalog], READ_ONLY, 0)

	# Find table row which matches obsmode

	call findmode (tp, obscol, obsmode, jrow)

	# Check for match, read and parse shape function if found

	if (jrow == 0) {
	    call printerr_str (badmode, obsmode)

	} else {
	    call syncolptr (tp, funccol, 2, cp)
	    call tbegtt (tp, cp, jrow, Memc[func], SZ_FNAME)

	    apshape = breakfunc (Memc[func])
	}

	call tbtclo (tp)
	call sfree (sp)
end

