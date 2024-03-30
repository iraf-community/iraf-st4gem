#* HISTORY *
#* B.Simon	25-May-94	Original
#* B.Simon	26-May-96	Use a single input spectrum
#* B.Simon	23-Jul-96	Renormalize after applying reddening
#* B.Simon	14-Aug-96	Add user selectable magnitude form
#* B.Simon	21-Nov-96	Add user selectable reddening law

# BUILDSPEC -- Build a spectral expression from its parts

procedure buildspec (spectrum, magnitude, magform, exptime, 
		     reddening, redlaw, newspec, maxch)

char	spectrum[ARB]	# i: user spectrum
char	magnitude[ARB]	# i: magnitude of synthetic spectrum
char	magform[ARB]	# i: magnitude units
real	exptime		# i: exposure time in seconds
real	reddening	# i: interstellar reddening E(B-V)
char	redlaw[ARB]	# i: name of reddening law
char	newspec[ARB]	# o: spectral expression replacing user spectrum
int	maxch		# i: maximum length of spectral expression
#--
int	ic
pointer	sp, redspec, factor, band
real	mag

string	nospectrum   "No spectrum specified"
string	badmagnitude "Parameter magnitude must start with a number"
string	badpassband  "Parameter magnitude does not contain passband"

bool	isblank(), is_simple()
int	ctor()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (redspec, SZ_FNAME, TY_CHAR)
	call salloc (factor, SZ_FNAME, TY_CHAR)
	call salloc (band, SZ_FNAME, TY_CHAR)

	if (isblank (spectrum))
	    call error (1, nospectrum)

	# Add reddening to spectrum

	if (reddening == 0.0) {
	    call strcpy (spectrum, Memc[redspec], SZ_FNAME)

	} else {
	    call sprintf (Memc[redspec], SZ_FNAME, "%s*ebmvx(%g,%s)")
	    call pargstr (spectrum)
	    call pargr (reddening)
	    call pargstr (redlaw)
	}

	# Renormalize spectral expression

	if (isblank (magnitude)) {
	    # Do not renormalize if magnitude is blank

	    call strcpy (Memc[redspec], newspec, maxch)

	} else {
	    # Renormalize spectrum to magnitude in passband

	    ic = 1
	    if (ctor (magnitude, ic, mag) <= 0)
		call error (1, badmagnitude)

	    call strcpy (magnitude[ic], Memc[factor], SZ_FNAME)
	    call strfix (Memc[factor])

	    if (isblank (Memc[factor]))
		call error (1, badpassband)

	    if (is_simple (Memc[factor])) {
		call sprintf (Memc[band], SZ_FNAME, "band(%s)")
		call pargstr (Memc[factor])
	    } else {
		call strcpy (Memc[factor], Memc[band], SZ_FNAME)
	    }

	    call sprintf (newspec, maxch, "rn(%s,%s,%g,%s)")
	    call pargstr (Memc[redspec])
	    call pargstr (Memc[band])
	    call pargr (mag)
	    call pargstr (magform)
	}

	# Multiply by exposure time

	if (exptime != 1.0) {
	    call sprintf (Memc[factor], SZ_FNAME, "*%g")
	    call pargr (exptime)

	    call strcat (Memc[factor], newspec, maxch)
	}

	call sfree (sp)
end
