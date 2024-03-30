include	<tbset.h>
include "simtwo.h"
include "object.h"
include	"otf.h"

define MAXCOL	5

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-May-95	Changed shape descriptor
#* B.Simon	02-Feb-98	fixed bug in copying flux to sflux

# RD_OBJECTS -- Read the list of objects from the input file

pointer procedure rd_objects (nwave, wave, thruput, input, spectrum, magband, 
			      magform, colnames, grftable, cmptable, skycoord,
			      exptime, apx, apy, psf, mw)

int	nwave		# i: length of wavelength and thruput arrays
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
char	input[ARB]	# i: object description file
char	spectrum[ARB]	# i: default spectrum
char	magband[ARB]	# i: magnitude passband
char	magform[ARB]	# i: magnitude form
char	colnames[ARB]	# i: input file column names (if table)
char	grftable[ARB]	# i: instrument graph table
char	cmptable[ARB]	# i: component name table
bool	skycoord	# i: use sky coordinates (ra & dec) for objects?
real	exptime		# i: exposure time
int	apx		# i: aperture x dimension
int	apy		# i: aperture y dimension
pointer	psf		# i: psf descriptor
pointer	mw		# i: coordinate transformation to aperture system
#--
bool	text
double	r[2], w[2], cd[2,2]
double	ra, dec, x, y
int	tp, nrow, ncol, irow, jrow, icol
pointer	cp[MAXCOL]
pointer	obj, sp, cname, units, spec, object, flux, sflux, ct
real	mag

string	raform    RA_UNITS
string	decform   DEC_UNITS
string	detform   DETECT_UNITS
string	objform   OBJ_UNITS

bool	isblank()
int	tbpsta(), word_find()
pointer	tbtopn(), mw_sctran(), breakfunc()
real	effstim()

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (cname, SZ_COLNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (object, SZ_LINE, TY_CHAR)
	call salloc (flux, nwave, TY_REAL)

	# Open object table, get column descriptors

	tp = tbtopn (input, READ_ONLY, NULL)

	text = tbpsta (tp, TBL_WHTYPE) == TBL_TYPE_TEXT
	nrow = tbpsta (tp, TBL_NROWS)
	ncol = tbpsta (tp, TBL_NCOLS)
	ncol = min (ncol, MAXCOL)

	Memc[cname] = EOS
	call amovki (NULL, cp, MAXCOL)
	do icol = 1, ncol {
	    if (! text) {
		if (word_find (icol, colnames, Memc[cname], SZ_COLNAME) == 0)
		    break
	    }

	    call syncolptr (tp, Memc[cname], icol, cp[icol])
	}

	# Allocate object structure

	call malloc (obj, SZ_OBJSTRUCT, TY_INT)
	call malloc (OBJ_XPOSARY(obj), nrow, TY_REAL)
	call malloc (OBJ_YPOSARY(obj), nrow, TY_REAL)
	call malloc (OBJ_FLUXARY(obj), nrow, TY_REAL)
	call malloc (OBJ_PWTARY(obj), nrow*OTF_NUMBER(psf), TY_REAL)
	call malloc (OBJ_SHPARY(obj), nrow, TY_INT)
	OBJ_NUMPSF(obj) = OTF_NUMBER(psf)

	# Read each line from input and compute object parameters

	sflux = NULL

	call mw_gwtermd (mw, r, w, cd, 2)
	ct = mw_sctran (mw, "world", "logical", 3)

	jrow = 0
	do irow = 1, nrow {
	    # Get object position

	    call tbegtd (tp, cp[1], irow, ra)
	    call tbegtd (tp, cp[2], irow, dec)

	    call tbcigt (cp[1], TBL_COL_UNITS, Memc[units], SZ_COLUNITS)
	    if (isblank (Memc[units])) {
		if (skycoord) {
		    call strcpy (raform, Memc[units], SZ_COLUNITS)
		} else {
		    call strcpy (detform, Memc[units], SZ_COLUNITS)
		}
	    }

	    call angtodeg (Memc[units], ra)

	    call tbcigt (cp[2], TBL_COL_UNITS, Memc[units], SZ_COLUNITS)
	    if (isblank (Memc[units])) {
		if (skycoord) {
		    call strcpy (decform, Memc[units], SZ_COLUNITS)
		} else {
		    call strcpy (detform, Memc[units], SZ_COLUNITS)
		}
	    }

	    call angtodeg (Memc[units], dec)

	    call ralimit (w[1], ra)
	    call mw_c2trand (ct, ra, dec, x, y)

	    # Do not include any object whose center is off the aperture

	    if (x < 0.5 || x > apx + 0.5 || y < 0.5 || y > apy + 0.5)
		next

	    jrow = jrow + 1
	    OBJ_XPOS(obj,jrow) = x
	    OBJ_YPOS(obj,jrow) = y

	    # Compute object flux

	    call tbegtr (tp, cp[3], irow, mag)

	    if (cp[4] == NULL) {
		Memc[spec] = EOS
	    } else {
		call tbegtt (tp, cp[4], irow, Memc[spec], SZ_LINE)
		call chkdefault (Memc[spec], SZ_LINE)
	    }

	    if (Memc[spec] != EOS) {
		call normspec (Memc[spec], exptime, mag, magband, magform, 
			      grftable, cmptable, nwave, wave, Memr[flux])

	    } else if (sflux != NULL) {
		call amovr (Memr[sflux], Memr[flux], nwave)

	    } else {
		call normspec (spectrum, exptime, mag, magband, magform, 
			      grftable, cmptable, nwave, wave, Memr[flux])

		call salloc (sflux, nwave, TY_REAL)
		call amovr (Memr[flux], Memr[sflux], nwave)
	    }


	    OBJ_FLUX(obj,jrow) = effstim (nwave, wave, thruput, 
					   Memr[flux], objform)

	    # Compute object psf weights

	    call psfweights (nwave, wave, thruput, Memr[flux], 
			     OTF_NUMBER(psf), OTF_WAVE(psf,1),
			     OBJ_PWEIGHT(obj,jrow,1))

	    # Compute object shape template

	    if (cp[5] == NULL) {
		Memc[object] = EOS
	    } else {
		call tbegtt (tp, cp[5], irow, Memc[object], SZ_LINE)
		call chkdefault (Memc[object], SZ_LINE)
	    }

	    if (Memc[object] == EOS) {
		OBJ_SHAPE(obj,jrow) = NULL
	    } else {
		OBJ_SHAPE(obj,jrow) = breakfunc (Memc[object])
	    }
	}

	OBJ_NUMBER(obj) = jrow

	call mw_ctfree (ct)
	call tbtclo (tp)
	call sfree (sp)
	return (obj)
end
