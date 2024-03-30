include	<tbset.h>

#* HISTORY *
#* B.Simon	14-Jun-94	original

# RDPHOTSTAT -- Get the statistical quantities frim the photmode table

procedure rdphotstat (tp, form, graftab, comptab, wave, nwave, 
		      stim, pivot, fwhm, nstat)

pointer	tp		# i: photometry table descriptor
char	form[ARB]	# i: output form
char	graftab[ARB]	# i: instrument graph table
char	comptab[ARB]	# i: component name table
real	wave[ARB]	# i: wavelength array 
int	nwave		# i: length of wavelength array
pointer	stim		# o: effective stimulus (must be freed by caller)
pointer	pivot		# o: pivot wavelength (must be freed by caller)
pointer	fwhm		# o: equivalent gaussian fwhm (must be freed by caller)
int	nstat		# o: length of output arrays
#--
int	nrow, irow, nmode, imode
pointer	sp, dt, fm, ob, oldmode, inform, mode, mode1, mode2
pointer	band1, band2, pcode, ocode1, ocode2

bool	streq()
int	tbpsta()
real	pivlam(), fwhmlam()

begin
	# Allocate memory for temporary strings and arrays

	call smark (sp)
	call salloc (oldmode, SZ_LINE, TY_CHAR)
	call salloc (inform, SZ_LINE, TY_CHAR)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (mode1, SZ_LINE, TY_CHAR)
	call salloc (mode2, SZ_LINE, TY_CHAR)
	call salloc (band1, nwave, TY_REAL)
	call salloc (band2, nwave, TY_REAL)
	call salloc (pcode, 2*SZ_LINE, TY_INT)

	ocode1 = pcode
	ocode2 = pcode + SZ_LINE

	# Allocate memory for output arrays

	nrow = tbpsta (tp, TBL_NROWS)

	call malloc (stim, 2*nrow, TY_REAL)
	call malloc (pivot, 2*nrow, TY_REAL)
	call malloc (fwhm, 2*nrow, TY_REAL)

	# Get pointers to data columns

	iferr {
	    call syncolptr (tp, "COUNTRATE", 1, dt)
	} then {
	    iferr {
		call syncolptr (tp, "DATUM", 1, dt)
	    } then {
		call syncolptr (tp, "FLUX", 1, dt)
	    }
	}

	call syncolptr (tp, "FORM", 2, fm)
	call syncolptr (tp, "OBSMODE", 3, ob)

	nstat = 0
	Memc[oldmode] = -1

	do irow = 1, nrow {
	    # Read data from table

	    call tbegtr (tp, dt, irow, Memr[stim+nstat])
	    call tbegtt (tp, fm, irow, Memc[inform], SZ_LINE)
	    call tbegtt (tp, ob, irow, Memc[mode], SZ_LINE)

	    if (streq (Memc[mode], Memc[oldmode])) {
		# If obsmode is not new, copy stats from last row

		call setphotform (Memc[inform], form, wave, nwave, 
				  Memr[band1], Memr[band2], nmode, 
				  Memr[stim+nstat], 1)

		if (imode > 1)
		    Memr[stim+nstat+1] = Memr[stim+nstat]

		do imode = 1, nmode {
		    Memr[pivot+nstat] = Memr[pivot+nstat-nmode]
		    Memr[fwhm+nstat] = Memr[fwhm+nstat-nmode]
		    nstat = nstat + 1
		}

	    } else {
		# Read obsmode, compute bandpass, and compute stats

		call strcpy (Memc[mode], Memc[oldmode], SZ_LINE)
		call splitexp (Memc[mode], nmode, Memc[mode1], 
			       Memc[mode2], SZ_LINE)

		call syncompile (Memc[mode1], Memi[ocode1], SZ_LINE)
		if (nmode > 1)
		    call syncompile (Memc[mode2], Memc[ocode2], SZ_LINE)

		call syncalc2 (graftab, comptab, Memc[mode1], Memc[mode2],
			       Memi[ocode1], Memi[ocode2], SZ_LINE, nwave,
			       wave, nmode, Memr[band1], Memr[band2])


		call setphotform (Memc[inform], form, wave, nwave, 
				  Memr[band1], Memr[band2], nmode, 
				  Memr[stim+nstat], 1)

		Memr[pivot+nstat] = pivlam (nwave, wave, Memr[band1])
		Memr[fwhm+nstat] = fwhmlam (nwave, wave, Memr[band1])
		nstat = nstat + 1

		if (nmode > 1) {
		    Memr[stim+nstat] = Memr[stim+nstat-1]
		    Memr[pivot+nstat] = pivlam (nwave, wave, Memr[band2])
		    Memr[fwhm+nstat] = fwhmlam (nwave, wave, Memr[band2])
		    nstat = nstat + 1
		}
	    }
	}

	call sfree (sp)
end
