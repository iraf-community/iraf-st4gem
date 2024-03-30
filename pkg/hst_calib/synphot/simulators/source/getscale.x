include	<tbset.h>
include	"simtwo.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original

# GETSCALE -- Retrieve the detector size and scale

procedure getscale (detcat, obsmode, apscale, apx, apy)

char	detcat[ARB]	# i: Catalog of detector information
char	obsmode[ARB]	# i: observation mode
double	apscale		# o: degrees per pixel
int	apx		# o: number pixels, x dimension
int	apy		# o: number pixels, y dimension
#--
int	jrow
pointer	sp, catalog, units, tp, cp

string	detform  APER_UNITS
string	obscol   "OBSMODE"
string	scalecol "SCALE"
string	xpixcol  "NX"
string	ypixcol  "NY"

string	badmode    "No match for obsmode found in detector catalog"

bool	isblank()
pointer	tbtopn()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catalog, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Open detector dimension  catalog

	call lastfile (detcat, Memc[catalog], SZ_FNAME)
	tp = tbtopn (Memc[catalog], READ_ONLY, 0)

	# Find table row which matches obsmode

	call findmode (tp, obscol, obsmode, jrow)

	# Check for match, copy values to output variables if matched

	if (jrow == 0) {
	    call printerr_str (badmode, obsmode)

	} else {
	    call syncolptr (tp, scalecol, 2, cp)
	    call tbegtd (tp, cp, jrow, apscale)

	    call tbcigt (cp, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	    # Use default units if units not found

	    if (isblank(Memc[units]))
		call strcpy (detform, Memc[units], SZ_COLUNITS)

	    call angtodeg (Memc[units], apscale)

	    call syncolptr (tp, xpixcol, 3, cp)
	    call tbegti (tp, cp, jrow, apx)
	    call syncolptr (tp, ypixcol, 4, cp)
	    call tbegti (tp, cp, jrow, apy)
	}


	call tbtclo (tp)
	call sfree (sp)
end
