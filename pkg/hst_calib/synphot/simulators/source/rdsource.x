include	<tbset.h>
include "simtwo.h"
include "source.h"
include	"otf.h"

define MAXCOL	5

# RD_SOURCE -- Read a list of spectral sources from 

procedure rd_source (nwave, wave, thruput, input, spectrum, magband, 
		     magform, colnames, grftable, cmptable, exptime, 
		     psf, lsf, mw, apshape, apscale, grating, dynrange, 
		     nsub, nx, ny, out)

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
real	exptime		# i: exposure time
pointer	psf		# i: psf descriptor
pointer	lsf		# i: lsf descriptor
pointer	mw		# i: coordinate transformation to aperture system
pointer	apshape		# i: aperture shape descriptor
double	apscale		# i: aperture scale
real	grating[ARB]	# i: grating dispersion parameters
real	dynrange	# i: dynamic range of object
int	nsub		# i: number of subpixels
int	nx		# i: first dimension of output buffer
int	ny		# i: second dimension of output buffer 
real	out[nx,ny]	# o: output buffer
#--
bool	text
double	r[2], w[2], cd[2,2]
double	ra, dec, x, y
int	tp, nrow, ncol, irow, icol, nweight, npx, npy
pointer	cp[MAXCOL]
pointer	src, sp, cname, units, spec, source, flux, weights, ct
real	mag

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
	call salloc (source, SZ_LINE, TY_CHAR)
	call salloc (flux, nwave, TY_REAL)

	if (lsf == NULL) {
	    nweight = OTF_NUMBER(psf)
	} else {
	    nweight = max (OTF_NUMBER(psf), OTF_NUMBER(lsf))
	}

	call salloc (weights, nweight, TY_REAL)

	# Open input table, get column descriptors

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

	# Read each line from input and compute source parameters

	call mw_gwtermd (mw, r, w, cd, 2)
	ct = mw_sctran (mw, "world", "logical", 3)

	do irow = 1, nrow {
	    # Allocate source structure

	    call malloc (src, SZ_SRCSTRUCT, TY_INT)
	    call malloc (SRC_SCALARS(src), SZ_SRCSCALARS, TY_REAL)
	    call malloc (SRC_WAVPTR(src), nwave, TY_REAL)
	    call malloc (SRC_SPECPTR(src), nwave, TY_REAL)

	    # Get source position

	    call tbegtd (tp, cp[1], irow, ra)
	    call tbegtd (tp, cp[2], irow, dec)

	    call tbcigt (cp[1], TBL_COL_UNITS, Memc[units], SZ_COLUNITS)
	    if (isblank (Memc[units]))
		call strcpy (detform, Memc[units], SZ_COLUNITS)

	    call angtodeg (Memc[units], ra)

	    call tbcigt (cp[2], TBL_COL_UNITS, Memc[units], SZ_COLUNITS)
	    if (isblank (Memc[units]))
		call strcpy (detform, Memc[units], SZ_COLUNITS)

	    call angtodeg (Memc[units], dec)

	    call ralimit (w[1], ra)
	    call mw_c2trand (ct, ra, dec, x, y)

	    SRC_XPOS(src) = x
	    SRC_YPOS(src) = y

	    # Compute source flux

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

	    } else {
		call normspec (spectrum, exptime, mag, magband, magform, 
			      grftable, cmptable, nwave, wave, Memr[flux])
	    }


	    SRC_FLUX(src) = effstim (nwave, wave, thruput, 
				     Memr[flux], objform)

	    SRC_NWAVE(src) = nwave
	    call amovr (wave, SRC_WAVE(src), nwave)
	    call amulr (thruput, Memr[flux], SRC_SPEC(src), nwave)

	    # Compute source point spread function

	    npx = OTF_NXPIX(psf)
	    npy = OTF_NYPIX(psf)
	    call malloc (SRC_PSFPTR(src), npx*npy, TY_REAL)

	    call psfweights (nwave, wave, thruput, SRC_SPEC(src),
			     OTF_NUMBER(psf), OTF_WAVE(psf,1),
			     Memr[weights])

	    SRC_NXPSF(src) = npx
	    SRC_NYPSF(src) = npy
	    call otfcalc (psf, x, y, nsub, OTF_NUMBER(psf), Memr[weights], 
			  npx, npy, SRC_PSF(src))

	    # Compute source line spread function

	    if (lsf == NULL) {
		SRC_NXLSF(src) = 0
		SRC_LSFPTR(src) = NULL

	    } else {
		npx = OTF_NXPIX(lsf)
		call malloc (SRC_LSFPTR(src), npx, TY_REAL)

		call psfweights (nwave, wave, thruput, SRC_SPEC(src),
				 OTF_NUMBER(lsf), OTF_WAVE(lsf,1),
				 Memr[weights])

		SRC_NXLSF(src) = npx
		call otfcalc (lsf, x, y, nsub, OTF_NUMBER(lsf), Memr[weights], 
			      npx, 1, SRC_LSF(src))

	    }

	    # Compute shape template

	    if (cp[5] == NULL) {
		Memc[source] = EOS
	    } else {
		call tbegtt (tp, cp[5], irow, Memc[source], SZ_LINE)
		call chkdefault (Memc[source], SZ_LINE)
	    }

	    if (Memc[source] == EOS) {
		SRC_SHAPE(src) = NULL
	    } else {
		SRC_SHAPE(src) = breakfunc (Memc[source])
	    }

	    # Add source to output array

	    call putsource (src, apshape, apscale, grating,
			    dynrange, nsub, nx, ny, out)

	    # Free spectral source structure

	    call free_func (SRC_SHAPE(src))
	    if (SRC_LSFPTR(src) != NULL)
		call mfree (SRC_LSFPTR(src), TY_REAL)

	    call mfree (SRC_PSFPTR(src), TY_REAL)
	    call mfree (SRC_SPECPTR(src), TY_REAL)
	    call mfree (SRC_WAVPTR(src), TY_REAL)
	    call mfree (SRC_SCALARS(src), TY_REAL)
	    call mfree (src, TY_INT)
	}

	call tbtclo (tp)
	call sfree (sp)

end
