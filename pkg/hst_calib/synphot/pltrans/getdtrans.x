include "dtrans.h"

#* HISTORY *
#* B.Simon	21-Jul-94	original

procedure getdtrans (spectrum, xmode, ymode, xform, yform, vzero, wavtab, 
		     graftab, comptab, dtrans)

char	spectrum[ARB]	# i: Target spectrum
char	xmode[ARB]	# u: x axis mode
char	ymode[ARB]	# u: y axis mode
char	xform[ARB]	# i: x axis form
char	yform[ARB]	# i: y axis form
char	vzero[ARB]	# i: Variable list
char	wavtab[ARB]	# i: Table containing wavelength array
char	graftab[ARB]	# i: Instrument graph table
char	comptab[ARB]	# i: Component name table
pointer	dtrans		# u: transformation data descriptor
#--
bool	nomode
int	nband, nxband, nyband, nwave, neff, istat, degree
real	v0
pointer	sp, xmode1, xmode2, ymode1, ymode2, spec, pcode, xcode1, xcode2
pointer	ycode1, ycode2, xband1, xband2, yband1, yband2, scode, sptr
pointer	wave, flux, xeff, yeff

string	badform  "Only the following forms are legal with no obsmode"
string	notspec  "Not a spectrum"

extern	getsynvar
int	is_count(), numlist(), numvzero(), nxtlist(), nxtvzero()
real	effstim2()
pointer	rdlist(), locpr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (xmode1, SZ_LINE, TY_CHAR)
	call salloc (xmode2, SZ_LINE, TY_CHAR)
	call salloc (ymode1, SZ_LINE, TY_CHAR)
	call salloc (ymode2, SZ_LINE, TY_CHAR)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (pcode, 5*SZ_LINE, TY_INT)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Replace empty obsmode with 1.0

	nomode = false
	if (xmode[1] == EOS) {
	    nomode = true
	    call strcpy ("1.0", xmode, SZ_FNAME)
	}

	if (ymode[1] == EOS) {
	    nomode = true
	    call strcpy ("1.0", ymode, SZ_FNAME)
	}

	if (nomode && (is_count (xform) == NO || is_count (yform) == NO))
	    call printerr_str (badform, "counts,obmag")

	# Compile the observation mode expressions

	nband = 0
	call splitexp (xmode, nxband, Memc[xmode1], Memc[xmode2], SZ_LINE)
	call splitexp (ymode, nyband, Memc[ymode1], Memc[ymode2], SZ_LINE)

	xcode1 = pcode + nband * SZ_LINE
	call syncompile (Memc[xmode1], Memi[xcode1], SZ_LINE)
	nband = nband + 1

	if (nxband == 1) {
	    xcode2 = xcode1
	} else {
	    xcode2 = pcode + nband * SZ_LINE
	    call syncompile (Memc[xmode2], Memi[xcode2], SZ_LINE)
	    nband = nband + 1
	}

	ycode1 = pcode + nband * SZ_LINE
	call syncompile (Memc[ymode1], Memi[ycode1], SZ_LINE)
	nband = nband + 1

	if (nyband == 1) {
	    ycode2 = ycode1
	} else {
	    ycode2 = pcode + nband * SZ_LINE
	    call syncompile (Memc[ymode2], Memi[ycode2], SZ_LINE)
	    nband = nband + 1
	}

	# Process the spectrum and vzero strings

	sptr = rdlist (spectrum)
	call rdvzero (vzero)

	# Process each spectrum

	wave = NULL
	while (nxtlist (sptr, Memc[spec], SZ_LINE) != EOF) {

	    scode = pcode + nband * SZ_LINE
	    call syncompile (Memc[spec], Memi[scode], SZ_LINE)

	    # Allocate wavealength set and descriptor 
	    # on first pass through loop

	    if (wave == NULL) {
		call getwavelen (wavtab, graftab, comptab, Memi[pcode], 
				 nband+1, SZ_LINE, wave, nwave)

		# Compute observation mode passbands

		call salloc (xband1, nwave, TY_REAL)
		call salloc (xband2, nwave, TY_REAL)
		call salloc (yband1, nwave, TY_REAL)
		call salloc (yband2, nwave, TY_REAL)
		call salloc (flux, nwave, TY_REAL)

		call syncalc2 (graftab, comptab, Memc[xmode1], Memc[xmode2], 
			       Memi[xcode1], Memi[xcode2], SZ_LINE, nwave, 
			       Memr[wave], nxband, Memr[xband1], Memr[xband2])

		call syncalc2 (graftab, comptab, Memc[ymode1], Memc[ymode2], 
			       Memi[ycode1], Memi[ycode2], SZ_LINE, nwave, 
			       Memr[wave], nyband, Memr[yband1], Memr[yband2])

		# Allocate descriptor

		neff = numlist (sptr) * numvzero ()

		if (dtrans == NULL) {
		    istat = 0

		    call malloc (dtrans, LEN_TRNSTRUCT, TY_STRUCT)
		    call malloc (xeff, neff, TY_REAL)
		    call malloc (yeff, neff, TY_REAL)

		} else {
		    istat = TRN_NEFF(dtrans)
		    neff = neff + istat

		    call realloc (xeff, neff, TY_REAL)
		    call realloc (yeff, neff, TY_REAL)
		}

		TRN_NEFF(dtrans) = neff
		TRN_XEFF(dtrans) = xeff
		TRN_YEFF(dtrans) = yeff
	    }

	    # Loop over each value of vzero

	    while (nxtvzero (v0) != EOF) {

		# Calculate effective stimulus 
		# of spectrum / passband combinations

		call syncalc (Memi[scode], SZ_LINE, locpr(getsynvar), 
			      nwave, Memr[wave], graftab, comptab, 
			      Memr[flux], degree)

		if (degree != 1)
		    call printerr_str (notspec, Memc[spec])

		Memr[xeff+istat] = effstim2 (xform, nwave, Memr[wave], nxband, 
					     Memr[xband1], Memr[xband2], 
					     Memr[flux])

		Memr[yeff+istat] = effstim2 (yform, nwave, Memr[wave], nyband, 
					     Memr[yband1], Memr[yband2], 
					     Memr[flux])

		istat = istat + 1
	    }
	}

	# Close files and release memory

	call clssyntab
	call freelist (sptr)

	if (wave != NULL)
	    call mfree (wave, TY_REAL)

	call sfree (sp)
end
