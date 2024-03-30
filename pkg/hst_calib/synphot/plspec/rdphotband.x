include	<tbset.h>

#* HISTORY *
#* B.Simon	14-Jun-94	original

# RDPHOTBAND -- Get the bandpasses from the obsmodes in a photometry table

procedure rdphotband (tp, graftab, comptab, wave, nwave, band, nband)

pointer	tp		# i: photometry table descriptor
char	graftab[ARB]	# i: instrument graph table
char	comptab[ARB]	# i: component name table
real	wave[ARB]	# i: wavelength array 
int	nwave		# i: length of wavelength array
pointer	band		# o: bandpass arrays (must be freed by caller)
int	nband		# o: number of bandpasses
#--
int	nrow, irow, nmode, iband
pointer	sp, cp, oldmode, mode, mode1, mode2, pcode, ocode1, ocode2

bool	streq()
int	tbpsta()

begin
	# Allocate memory for temporary strings and arrays

	call smark (sp)
	call salloc (oldmode, SZ_LINE, TY_CHAR)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (mode1, SZ_LINE, TY_CHAR)
	call salloc (mode2, SZ_LINE, TY_CHAR)
	call salloc (pcode, 2*SZ_LINE, TY_INT)

	ocode1 = pcode
	ocode2 = pcode + SZ_LINE

	# Count number of passbands in table

	nband = 0
	Memc[oldmode] = EOS
	nrow = tbpsta (tp, TBL_NROWS)
	call syncolptr (tp, "OBSMODE", 3, cp)

	do irow = 1, nrow {
	    call tbegtt (tp, cp, irow, Memc[mode], SZ_LINE)
	    if (streq (Memc[mode], Memc[oldmode]))
		next

	    call strcpy (Memc[mode], Memc[oldmode], SZ_LINE)
	    call splitexp (Memc[mode], nmode, Memc[mode1], 
			   Memc[mode2], SZ_LINE)

	    nband = nband + nmode
	}

	# Read obsmodes and convert into passbands

	iband = 0
	Memc[oldmode] = EOS
	call malloc (band, nwave*(nband+1), TY_REAL)

	do irow = 1, nrow {
	    call tbegtt (tp, cp, irow, Memc[mode], SZ_LINE)
	    if (streq (Memc[mode], Memc[oldmode]))
		next

	    call strcpy (Memc[mode], Memc[oldmode], SZ_LINE)
	    call splitexp (Memc[mode], nmode, Memc[mode1], 
			   Memc[mode2], SZ_LINE)

	    call syncompile (Memc[mode1], Memi[ocode1], SZ_LINE)
	    if (nmode > 1)
		call syncompile (Memc[mode2], Memc[ocode2], SZ_LINE)

	    call syncalc2 (graftab, comptab, Memc[mode1], Memc[mode2],
			   Memi[ocode1], Memi[ocode2], SZ_LINE, nwave,
			   wave, nmode, Memr[band+iband*nwave], 
			   Memr[band+(iband+1)*nwave])

	    iband = iband + nmode
	}

	call sfree (sp)
end
