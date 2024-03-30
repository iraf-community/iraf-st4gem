include	<tbset.h>

#* HISTORY *
#* B.Simon	17-Feb-95	original

# OTFOPEN -- Open the otf images associated with the observation mode

procedure otfopen (otfcat, obsmode, im, waveotf, numotf)

char	otfcat[ARB]	# i: catalog of otf file names
char	obsmode[ARB]	# i: observation mode
pointer	im		# o: point spread function image descriptors
pointer	waveotf		# o: associated wavelengths for otfs
int	numotf		# o: number of point spred functions
#--
int	maxkey, nkey, nrow, irow, ic, idx
pointer	sp, catalog, catmode, key, file, tp, ob, wv, ps
real	wavelen

string	obscol  "OBSMODE"
string	wavecol "WAVELENGTH"
string	otfcol  "FILENAME"
string	nofile   "PSF or LSF catalog name is blank. Check simcatp"
string	notfound "PSF or LSF catalog not found"
string	badmode  "No match for obsmode found in PSF or LSF catalog"

bool	isblank()
int	access(), imaccess(), tbpsta(), word_fetch(), word_match()
pointer	tbtopn(), immap()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catalog, SZ_FNAME, TY_CHAR)
	call salloc (catmode, SZ_FNAME, TY_CHAR)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (file, SZ_FNAME, TY_CHAR)

	# Check for missing catalog

	if (isblank (otfcat))
	    call printerr_str (nofile, otfcat)

	# Open otf catalog

	call lastfile (otfcat, Memc[catalog], SZ_FNAME)
	if (access (Memc[catalog], 0, 0) == NO)
	    call printerr_str (notfound, Memc[catalog])

	# If the "catalog file" is a single image, open that
	# image as the sole otf. Otherwise, open the catalog and
	# search for the match otf's

	if (imaccess (Memc[catalog], READ_ONLY) == YES &&
	    imaccess (Memc[catalog], NEW_FILE) == YES    ) {

	    numotf = 1
	    call malloc (im, numotf, TY_INT)
	    call malloc (waveotf, numotf, TY_REAL)

	    Memi[im] = immap (Memc[catalog], READ_ONLY, 0)
	    Memr[waveotf] = 5500.	# dummy wavelength

	} else {
	    tp = tbtopn (Memc[catalog], READ_ONLY, 0)
	    nrow = tbpsta (tp, TBL_NROWS)

	    call syncolptr (tp, obscol, 1, ob)
	    call syncolptr (tp, wavecol, 2, wv)
	    call syncolptr (tp, otfcol, 3, ps)

	    # Allocate memory to hold otf descriptors

	    call malloc (im, nrow, TY_INT)
	    call malloc (waveotf, nrow, TY_REAL)

	    # Search catalog for rows matching obsmode

	    numotf = 0
	    maxkey = 0
	    do irow = 1, nrow {
		call tbegtt (tp, ob, irow, Memc[catmode], SZ_FNAME)

		# Check each keyword in catalog obsmode against input obsmode

		ic = 1
		nkey = 0
		while (word_fetch (Memc[catmode], ic, Memc[key], SZ_FNAME) > 0) {
		    if (word_match (Memc[key], obsmode) == 0) {
			nkey = 0
			break
		    }

		    nkey = nkey + 1
		}

		# Save otf information if match

		if (nkey >= maxkey && nkey > 0) {
		    # Reinitialize list if better match found

		    if (nkey > maxkey) {
			do idx = 0, numotf-1
			    call  imunmap (Memi[im+idx])

			numotf = 0
			maxkey = nkey
		    }

		    # Copy wavelength and image to output arrays

		    call tbegtt (tp, ps, irow, Memc[file], SZ_FNAME)
		    call tbegtr (tp, wv, irow, wavelen)

		    # Insert in sorted order of wavelengths

		    for (idx = numotf; idx > 0; idx = idx - 1) {
			if (wavelen >= Memr[waveotf+idx-1])
			    break

			Memi[im+idx] = Memi[im+idx-1]
			Memr[waveotf+idx] = Memr[waveotf+idx-1]
		    }

		    Memi[im+idx] = immap (Memc[file], READ_ONLY, 0)
		    Memr[waveotf+idx] = wavelen
		    numotf = numotf + 1
		}
	    }

	    # Release memory and close catalog

	    call tbtclo (tp)

	    if (numotf == 0) {
		call mfree (im, TY_INT)
		call mfree (waveotf, TY_REAL)
		call printerr_str (badmode, obsmode)

	    } else {
		call realloc (im, numotf, TY_INT)
		call realloc (waveotf, numotf, TY_REAL)
	    }
	}

	call sfree (sp)
end
