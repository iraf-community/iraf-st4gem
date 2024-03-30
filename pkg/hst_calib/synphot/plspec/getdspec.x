include "dspec.h"

#* HISTORY *
#* B.Simon	08-Jun-94	original

# GETDSPEC -- Get spectra to plot and save in descriptor

procedure getdspec (obsmode, spectrum, form, vzero, 
		    wavtab, graftab, comptab, dspec)

char	obsmode[ARB]	# i: Instrument observation mode
char	spectrum[ARB]	# i: Target spectrum
char	form[ARB]	# i: Form of output spectrum
char	vzero[ARB]	# i: Variable list
char	wavtab[ARB]	# i: Table containing wavelength array
char	graftab[ARB]	# i: Instrument graph table
char	comptab[ARB]	# i: Component name table
pointer	dspec		# o: spectrum descriptor
#--
int	nband, nwave, nflux, nstat, degree, done, count
real	v0
pointer	sp, mode, mode1, mode2, spec, pcode, scode, ocode1, ocode2
pointer	optr, sptr, wave, band1, band2, flux, stim, pivot, fwhm

string	notspec  "Not a spectrum"

extern	getsynvar
int	is_count(), numlist(), numvzero(), nxtlist(), nxtvzero(), phottoany()
real	effstim2(), pivlam2(), fwhmlam2()
pointer	rdlist(), locpr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (mode1, SZ_LINE, TY_CHAR)
	call salloc (mode2, SZ_LINE, TY_CHAR)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (pcode, 3*SZ_LINE, TY_INT)

	scode = pcode
	ocode1 = pcode + SZ_LINE
	ocode2 = pcode + 2 * SZ_LINE

	# Replace empty obsmode with 1.0

	if (obsmode[1] == EOS) {
	    call strcpy ("1.0", Memc[mode], SZ_LINE)
	} else {
	    call strcpy (obsmode, Memc[mode], SZ_LINE)
	}

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Process the obsmode, spectrum,  and vzero strings

	optr = rdlist (Memc[mode])
	sptr = rdlist (spectrum)
	call rdvzero (vzero)

	# Loop over each expression

	dspec = NULL
	count = is_count (form)

	while (nxtlist (optr, Memc[mode], SZ_LINE) != EOF) {
	    call splitexp (Memc[mode], nband, Memc[mode1], 
			   Memc[mode2], SZ_LINE)

	    call syncompile (Memc[mode1], Memi[ocode1], SZ_LINE)
	    if (nband == 1) {
	        Memi[ocode2] = -1
	    } else {
		call syncompile (Memc[mode2], Memi[ocode2], SZ_LINE)
	    }

	    # Compute the obsmode bandpass (except the first)

	    if (dspec != NULL) {
		call syncalc2 (graftab, comptab, Memc[mode1], Memc[mode2], 
			       Memi[ocode1], Memi[ocode2], SZ_LINE, nwave, 
			       Memr[wave], nband, Memr[band1], Memr[band2])
	    }

	    # Process each target / spectrum

	    while (nxtlist (sptr, Memc[spec], SZ_LINE) != EOF) {
		call syncompile (Memc[spec], Memi[scode], SZ_LINE)

		# The first obsmode is processed here since we have to wait
		# until first spectrum is read to compute the wavelength set

		# Allocate descriptor on first pass through loop

		if (dspec == NULL) {
		    call getwavelen (wavtab, graftab, comptab, Memi[pcode], 
				     nband+1, SZ_LINE, wave, nwave)

		    nflux = numlist(optr) * numlist (sptr) * numvzero ()
		    nstat = 0

		    call salloc (band1, nwave, TY_REAL)
		    call salloc (band2, nwave, TY_REAL)

		    call malloc (dspec, LEN_SPCSTRUCT, TY_STRUCT)
		    call malloc (flux, nwave*nflux, TY_REAL)
		    call malloc (stim, nflux, TY_REAL)
		    call malloc (pivot, nflux, TY_REAL)
		    call malloc (fwhm, nflux, TY_REAL)


		    SPC_NWAVE(dspec) = nwave
		    SPC_NFLUX(dspec) = nflux
		    SPC_WAVE(dspec) = wave
		    SPC_FLUX(dspec) = flux
		    SPC_STIM(dspec) = stim
		    SPC_PIVOT(dspec) = pivot
		    SPC_FWHM(dspec) = fwhm

		    call syncalc2 (graftab, comptab, Memc[mode1], Memc[mode2], 
				   Memi[ocode1], Memi[ocode2], SZ_LINE, nwave, 
				   Memr[wave], nband, Memr[band1], Memr[band2])
		}

		# Loop over each value of vzero

		while (nxtvzero (v0) != EOF) {

		    # Calculate flux of spectrum / passband combination

		    call syncalc (Memi[scode], SZ_LINE, locpr(getsynvar), 
				  nwave, Memr[wave], graftab, comptab, 
				  Memr[flux], degree)

		    if (degree != 1)
			call printerr_str (notspec, Memc[spec])

		    # Compute statistics

		    if (obsmode[1] == EOS && count == NO) {
			Memr[stim+nstat] = INDEFR

		    } else {
			Memr[stim+nstat] = effstim2 (form, nwave, Memr[wave], 
						     nband, Memr[band1], 
						     Memr[band2], Memr[flux])
		    }

		    Memr[pivot+nstat] = pivlam2 (nwave, Memr[wave], nband,
						 Memr[band1], Memr[band2])

		    Memr[fwhm+nstat] = fwhmlam2 (nwave, Memr[wave], nband,
						 Memr[band1], Memr[band2])

		    nstat = nstat + 1

		    # Multiply flux by bandpass to get total flux 

		    call bandmul (nwave, nband, Memr[band1], Memr[band2], 
				  Memr[flux])

		    # Convert flux to output units

		    if (nband == 1) {
			done = phottoany (form, nwave, Memr[wave], Memr[flux])
		    }

		    flux = flux + nwave
		}
	    }
	}

	# Close files and release memory

	call clssyntab
	call freelist (optr)
	call freelist (sptr)

	call sfree (sp)
end
