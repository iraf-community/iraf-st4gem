include <tbset.h>

# GRIDFUNC -- Interpolate in a grid of spectra

procedure gridfunc (list, place, nwave, wave, spec)

char	list[ARB]	# i: list of filenames
real	place		# i: place to interpolate at
int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	spec[ARB]	# o: output spectrum
#--
int	ilo, ihi, nrow
pointer	sp, tp, cp, loname, hiname, lospec, hispec
real	loval, hival

string	noextrap  "Cannot extrapolate from list of spectra"

int	tbpsta(), opnsyntab()
errchk	opnsyntab, syncolptr, tbeggt, rdspec, synphoterr

begin
	# Allocate memory for temporary strings and arrays

	call smark (sp)
	call salloc (loname, SZ_FNAME, TY_CHAR)
	call salloc (hiname, SZ_FNAME, TY_CHAR)
	call salloc (lospec, nwave, TY_REAL)
	call salloc (hispec, nwave, TY_REAL)

	# Check arguments

	tp = opnsyntab (list)
	nrow = tbpsta (tp, TBL_NROWS)

	ilo = int (place)
	ihi = min (ilo + 1, nrow)

	if (ilo < 1 || ilo > nrow)
	    call synphoterr (noextrap, "grid")

	# Read spectra which bound interpolating point

	call syncolptr (tp, "FILENAME", 1, cp)
	call tbegtt (tp, cp, ilo, Memc[loname], SZ_FNAME)
	call tbegtt (tp, cp, ihi, Memc[hiname], SZ_FNAME)

	call rdspec (Memc[loname], nwave, wave, Memr[lospec])
	call rdspec (Memc[hiname], nwave, wave, Memr[hispec])

	# Interpolate to produce output spectrum

	loval = ilo
	hival = ihi

	call funinterp (loval, hival, place, nwave, 
			Memr[lospec], Memr[hispec], spec)

	call sfree (sp)
end
