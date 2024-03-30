#* HISTORY *
#* B.Simon	14-Jun-94	original

# RDPHOTWAVE -- Extract a wavelength array from the photometry table

procedure rdphotwave (tp, wavtab, graftab, comptab, wave, nwave)

pointer	tp		# i: photometry table descriptor
char	wavtab[ARB]	# i: table containing wavelength array
char	graftab[ARB]	# i: instrument graph table
char	comptab[ARB]	# i: component name table
pointer	wave		# o: wavelength array (must be freed by caller)
int	nwave		# o: length of wavelength array
#--
int	nband
pointer	sp, mode, mode1, mode2, pcode, ocode1, ocode2, cp

begin
	# Allocate memory for tempo
	call smark (sp)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (mode1, SZ_LINE, TY_CHAR)
	call salloc (mode2, SZ_LINE, TY_CHAR)
	call salloc (pcode, 2*SZ_LINE, TY_INT)

	ocode1 = pcode
	ocode2 = pcode + SZ_LINE

	# Read observation mode from first row of table

	call syncolptr (tp, "OBSMODE", 0, cp)
	call tbegtt (tp, cp, 1, Memc[mode], SZ_LINE)

	# Derive a wavelength array from it

	call splitexp (Memc[mode], nband, Memc[mode1], Memc[mode2], SZ_LINE)

	call syncompile (Memc[mode1], Memi[ocode1], SZ_LINE)
	if (nband > 1)
	    call syncompile (Memc[mode2], Memc[ocode2], SZ_LINE)

	call getwavelen (wavtab, graftab, comptab, Memi[pcode],
			 nband, SZ_LINE, wave, nwave)

	call sfree (sp)
end
