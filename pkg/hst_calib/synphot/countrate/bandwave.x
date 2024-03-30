include	<tbset.h>

#* HISTORY *
#* B.Simon	25-Mar-93	original
#* B.Simon	26-May-94	modified to call rdwave and return flag
#* B.Simon	02-Jun-94	modified to strip band() from mode string
#* B.Simon	20-Mar-95	added call to modefile()
#* B.Simon	26-Mar-96	added call to calcwave()
#* B.Simon	20-Aug-96	pass wavetab instead of wavecat

# BANDWAVE -- Compute the proper wavelength set for a given bandpass

bool procedure bandwave (wavetab, band, cenwave, nwave, wave)

char	wavetab[ARB]	# i: wavelength table name
char	band[ARB]	# i: instrument mode string
real	cenwave		# i: central wavelength
int	nwave		# o: number of wavelengths in set
pointer	wave		# o: pointer to array of wavelengths
#--
bool	status
int	ic, index, iwave, first, last, maxwave[11]
pointer	sp, tp, key, list, owave

data	maxwave	/ 0, 0, 0, 0, 0, 0, 0, 500, 0, 1024, 0 /

string	nowavetab  "Warning: cannot read wavelength table, using default\n"
string	noinst     "bandwave: no instrument found in instrument mode string"
string	instrume   "fgs,acs,wfp2,hsp,wfpc,foc,fos,hrs,nicmos,stis,wfc3"

int	findwave(), word_fetch(), word_match(), tbpsta()
pointer	tbtopn()

begin
	# Allocate dynamic memory for temporary strings

	call smark (sp)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (list, SZ_FNAME, TY_CHAR)

	# Set procedure return status

	status = true

	# Extract list of obsmode keywords from mode string

	call getnaked (band, Memc[list], SZ_FNAME)

	# Read wavelength table

	if (wavetab[1] == '(') {
	    call calcwave (wavetab, wave, nwave)

	} else ifnoerr (tp = tbtopn (wavetab, READ_ONLY, NULL)) {
	    nwave = tbpsta (tp, TBL_NROWS)
	    call malloc (wave, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)

	} else {
	    call eprintf (nowavetab)
	    status = false
	    wave = NULL
	    nwave = 0
	}

	if (status && ! IS_INDEFR(cenwave)) {
	    # Extract region containing central wavelength

	    ic = 1
	    index = 0
	    while (word_fetch (Memc[list], ic, Memc[key], SZ_FNAME) > 0) {
		index = word_match (Memc[key], instrume)
		if (index > 0)
		    break
	    }

	    if (index == 0)
		call error (1, noinst)

	    if (maxwave[index] > 0) {
		iwave = findwave (cenwave, Memr[wave], nwave)
		first = max (1, iwave - maxwave[index] / 2)
		last = min (nwave, first + maxwave[index] - 1)

		owave = wave
		nwave = last - first + 1

		call malloc (wave, nwave, TY_REAL)
		call amovr (Memr[owave+first-1], Memr[wave], nwave)
		call mfree (owave, TY_REAL)
	    }
	}

	call sfree (sp)
	return (status)

end
