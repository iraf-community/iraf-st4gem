include "dsphot.h"

#* HISTORY *
#* B.Simon	10-Jun-94	original

# GETDSPHOT -- Get spectrophotemetry files to plot and save in descriptor

procedure getdsphot (spfile, form, wavtab, dsphot)

char	spfile[ARB]	# i: Spectrophotometry file
char	form[ARB]	# i: Form of spectrum
char	wavtab[ARB]	# i: Wavelength table name
pointer	dsphot		# o: specrophotometry descriptor
#--
int	nflux, nwave
pointer	sp, file, sptr, wave, flux, err, fwhm

string	fluxcol  "FLUX"

int	numlist(), nxtlist()
pointer	rdlist()

begin
	call smark (sp)
	call salloc (file, SZ_LINE, TY_CHAR)

	dsphot = NULL
	call inisyntab

	sptr = rdlist (spfile)
	nflux = numlist (sptr)

	call listwave (wavtab, sptr, fluxcol, wave, nwave)

	while (nxtlist (sptr, Memc[file], SZ_LINE) != EOF) {
	    if (dsphot == NULL) {
		call malloc (flux, nflux*nwave, TY_REAL)
		call malloc (err, nflux*nwave, TY_REAL)
		call malloc (fwhm, nflux*nwave, TY_REAL)

		call malloc (dsphot, LEN_SPTSTRUCT, TY_STRUCT)
		SPT_NWAVE(dsphot) = nwave
		SPT_NFLUX(dsphot) = nflux
		SPT_WAVE(dsphot) = wave
		SPT_FLUX(dsphot) = flux
		SPT_ERROR(dsphot) = err
		SPT_FWHM(dsphot) = fwhm
	    }

	    call rdsphot (Memc[file], nwave, Memr[wave], 
			  Memr[flux], Memr[err], Memr[fwhm])

	    call errform (form, nwave, Memr[wave], Memr[flux], Memr[err])

	    flux = flux + nwave
	    err = err + nwave
	    fwhm = fwhm + nwave
	}

	call clssyntab
	call freelist (sptr)
	call sfree (sp)
end
