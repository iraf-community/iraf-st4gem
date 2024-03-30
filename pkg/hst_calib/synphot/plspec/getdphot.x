include	"dphot.h"

#* HISTORY *
#* B.Simon	14-Jun-94	original

# GETDPHOT -- Get photometry table info  to plot and save in descriptor

procedure getdphot (pfile, form, wavtab, graftab, comptab, dphot)

char	pfile[ARB]	# i: photometry file
char	form[ARB]	# i: output form
char	wavtab[ARB]	# i: Table containing wavelength array
char	graftab[ARB]	# i: Instrument graph table
char	comptab[ARB]	# i: Component name table
pointer	dphot		# o: photometry descriptor
#--
int	nstat, nband, maxstat, maxband, nwave, iband, istat
pointer	sp, phottab, sptr, wave, band, stim, pivot, fwhm
pointer	tp, pband, pstim, ppiv, pfwhm

string	badphot  "Could not open photometry table"

int	nxtlist()
pointer	rdlist(), tbtopn()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (phottab, SZ_LINE, TY_CHAR)

	call inisyntab

	sptr = rdlist (pfile)

	# Initialize array pointers and length of arrays

	wave = NULL
	band = NULL
	stim = NULL
	pivot = NULL
	fwhm = NULL

	nstat = 0
	nband = 0
	maxstat = 0
	maxband = 0

	while (nxtlist (sptr, Memc[phottab], SZ_LINE) != EOF) {
	    # Open photometry table

	    iferr {
		tp = tbtopn (Memc[phottab], READ_ONLY, NULL)
	    } then {
		call synphotwarn (badphot, Memc[phottab])
		next
	    }

	    # Get wavelength set from first photometry table

	    if (wave == NULL)
		call rdphotwave (tp, wavtab, graftab, 
				 comptab, wave, nwave)

	    # Read the bandpasses from the table and store in array

	    call rdphotband (tp, graftab, comptab, Memr[wave], 
			     nwave, pband, iband)

	    maxband = maxband + iband * nwave
	    call realloc (band, maxband, TY_REAL)

	    call amovr (Memr[pband], Memr[band+nband*nwave], iband*nwave)
	    nband = nband + iband

	    call mfree (pband, TY_REAL)

	    # Read the statistics from the table and store in array

	    call rdphotstat (tp, form, graftab, comptab, Memr[wave], nwave, 
			     pstim, ppiv, pfwhm, istat)

	    maxstat = maxstat + istat
	    
	    call realloc (stim, maxstat, TY_REAL)
	    call realloc (pivot, maxstat, TY_REAL)
	    call realloc (fwhm, maxstat, TY_REAL)

	    call amovr (Memr[pstim], Memr[stim+nstat], istat)
	    call amovr (Memr[ppiv], Memr[pivot+nstat], istat)
	    call amovr (Memr[pfwhm], Memr[fwhm+nstat], istat)

	    nstat = nstat + istat

	    call mfree (pstim, TY_REAL)
	    call mfree (ppiv, TY_REAL)
	    call mfree (pfwhm, TY_REAL)

	    call tbtclo (tp)
	}

	# Save arrays in photometry descriptor

	if (wave == NULL) {
	    dphot = NULL

	} else {
	    call malloc (dphot, LEN_PHTSTRUCT, TY_STRUCT)
	    PHT_NWAVE(dphot) = nwave
	    PHT_NBAND(dphot) = nband
	    PHT_WAVE(dphot) = wave
	    PHT_BAND(dphot) = band
	    PHT_NSTAT(dphot) = nstat
	    PHT_STIM(dphot) = stim
	    PHT_PIVOT(dphot) = pivot
	    PHT_FWHM(dphot) = fwhm
	}

	call clssyntab
	call freelist (sptr)
	call sfree (sp)
end
