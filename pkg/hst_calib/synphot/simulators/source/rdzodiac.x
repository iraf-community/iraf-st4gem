include	<tbset.h>

#* HISTORY *
#* B.Simon	14-Apr-95	original

# RDZODIAC -- Look up flux of zodiacal light in table

procedure rdzodiac (zodtab, ra, dec, jd, zflux)

char	zodtab[ARB]	# i: name of zodiacal light table
double	ra		# i: right ascension of target
double	dec		# i: declination of target
double	jd		# i: julian date of observation
real	zflux		# o: zodiacal flux in 10th mag stars / deg ** 2
#--
double	lh, lt, bt, lz, bz
int	ncol, nrow, icol, irow, ilo, ihi
pointer	tp, cp, cplo, cphi
real	lzlo, lzhi, bzlo, bzhi, dl, db, al, bl, ab, bb
real	z1, z2, z3, z4

string	nobrak  "Cannot bracket target position in table"

int	tbpsta()
pointer	tbtopn(), tbcnum()

begin
	# Compute ecliptic longitude and latitude of sun and target

	call sunlong (jd, lh)
	call ecliptic (ra, dec, lt, bt)

	# Compute angles between sun and target

	call ralimit (lh, lt)
	lt = abs (lt - lh)
	bt = abs (bt)

	# Find columns and rows which bracket these angles

	tp = tbtopn (zodtab, READ_ONLY, NULL)

	ncol = tbpsta (tp, TBL_NCOLS)
	nrow = tbpsta (tp, TBL_NROWS)

	# Find latitude bracket
	# First row in table contains latitudes at which 
	# flux values are tabulated

	cplo = NULL
	cphi = NULL

	bzlo = -1.0
	bzhi = 91.0

	do icol = 2, ncol {
	    cp = tbcnum (tp, icol)
	    call tbegtd (tp, cp, 1, bz)

	    if (bz <= bt) {
		if (bz > bzlo) {
		    bzlo = bz
		    cplo = cp
		}
	    }

	    if (bz >= bt) {
		if (bz < bzhi) {
		    bzhi = bz
		    cphi = cp
		}
	    }
	}

	if (cplo == NULL || cphi == NULL)
	    call printerr_str (nobrak, zodtab)

	# Find longitude bracket
	# First column of table contains longitudes at which
	# flux values are tabulated

	ilo = 0
	ihi = 0

	lzlo = -1.0
	lzhi = 181.0

	cp = tbcnum (tp, 1)

	do irow = 2, nrow {
	    call tbegtd (tp, cp, irow, lz)

	    if (lz <= lt) {
		if (lz > lzlo) {
		    lzlo = lz
		    ilo = irow
		}
	    }

	    if (lz >= lt) {
		if (lz < lzhi) {
		    lzhi = lz
		    ihi = irow
		}
	    }
	}

	if (ilo == 0 || ihi == 0)
	    call printerr_str (nobrak, zodtab)

	# Compute interpolation factors:
	# al, bl = longitude factors
	# ab, bb = latitude factors

	dl = lzhi - lzlo

	if (dl == 0.0) {
	    bl = 1.0
	} else {
	    bl = (lt - lzlo) / dl
	}

	al = 1.0 - bl

	db = bzhi - bzlo
	if (db == 0.0) {
	    bb = 1.0
	} else {
	    bb = (bt - bzlo) / db
	}

	ab = 1.0 - bb

	# Interpolate in table for flux value

	call tbegtr (tp, cplo, ilo, z1)
	call tbegtr (tp, cplo, ihi, z2)
	call tbegtr (tp, cphi, ilo, z3)
	call tbegtr (tp, cphi, ihi, z4)

	zflux = al * ab * z1 + bl * ab * z2 + al * bb * z3 + bl * bb * z4

	call tbtclo (tp)
end
