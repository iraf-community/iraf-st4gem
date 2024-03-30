include	"mac.h"

# RENORM -- Renormalize the spectrum
#
# Sep 1993 Howard Bushouse: For form=counts or obmag, call RATE, not EFFSTIM.
# Sep 1993 HB: For form=counts or obmag, descale requested normalization
#              value by HST area.

procedure renorm( script, grtbl, cmptbl, form, nwv, wv, iw, spec)

char	script[ARB]	# Command script
char 	grtbl[ARB]	# i: graph table name
char 	cmptbl[ARB]	# i: component table name
char 	form[ARB]	# i: form (units) of input spectrum
int 	nwv		# i: number of wavelengths
pointer	wv		# i: pointer to array of wavelengths
int	iw		# io: position in script
real 	spec[ARB]	# io: renormalized spectrum

int 	ip, nchar
int 	ctor(), ctowrd(), strsearch()
pointer	filt
real 	want, have, add, factor, hstarea
real 	effstim(), rate(), clgetr()
char	rnform[SZ_FNAME], mode[SZ_LINE], errmsg[SZ_LINE], rnval[SZ_FNAME]
bool 	streq(), status

string 	badrenorm 	"Invalid Renormalization."

begin

	# Mode
	nchar = ctowrd( script, iw, mode, SZ_LINE)

	# Form
	nchar = ctowrd( script, iw, rnform, SZ_FNAME)

	# Value
	nchar = ctowrd( script, iw, rnval, SZ_FNAME)
	ip = 1
	nchar = ctor( rnval, ip, want)

	# If rnform=counts or obmag, descale requested value by HST area
	if (streq(rnform, "counts")  || streq(rnform, "COUNTS") ||
	    streq(rnform, "photpix") || streq(rnform, "PHOTPIX") ) {
	    hstarea = clgetr("area")
	    want = want / hstarea
	} else if (streq(rnform, "obmag") || streq(rnform, "OBMAG") ) {
	    hstarea = clgetr("area")
	    want = want + 2.5 * alog10(hstarea)
	}

	ip = 1
	call compband( mode, ip, grtbl, cmptbl, nwv, wv, filt)

	call specform( nwv, Memr[wv], spec, form, spec, rnform, status )

	# count rate
	if (streq(rnform, "counts")  || streq(rnform, "COUNTS")  ||
	    streq(rnform, "photpix") || streq(rnform, "PHOTPIX") ||
	    streq(rnform, "obmag")   || streq(rnform, "OBMAG") ) {
	    have = rate( nwv, Memr[wv], Memr[filt], spec, rnform)
	    if (streq(rnform, "obmag") || streq(rnform, "OBMAG") )
		have = -2.5 * alog10(have)
	} else
	    have = effstim( nwv, Memr[wv], Memr[filt], spec, rnform )

	if ( strsearch( rnform, "mag") > 0 ) {
	   add =  want - have
	   call baddkr(spec, add, spec, nwv)

	} else if( have > 0. ) {
	   factor = want / have
	   call bmulkr( spec, factor, spec, nwv )

	} else {
	   call sprintf( errmsg, SZ_LINE, badrenorm)
	   call error( 1, errmsg)
	}
	call specform( nwv, Memr[wv], spec, rnform, spec, form, status )

	# Free allocated memory
	call mfree( filt, TY_REAL)
end
