include	<tbset.h>

define	SZ_WAVE		1200
define	SZ_UNIT		10

# MAKEWAVE -- Read a wavelength set from a table or create a dummy set

procedure makewave (wavetable, nwave, wave, units)

char	wavetable[ARB]	# i: name of wavelenght table
int	nwave		# o: number of wavelengths
pointer	wave		# o: wavelength array
pointer	units		# o: wavelength units
#--
pointer	sp, tp, table

bool	streq()
int	tbpsta()
pointer	tbtopn()

begin
	# Lower case table name to and remove blanks to simplify check

	call smark (sp)
	call salloc (table, SZ_FNAME, TY_CHAR)

	call strcpy (wavetable, Memc[table], SZ_FNAME)
	call strfix (Memc[table])

	# If no table name, create dummy table. The first value in the
	# dummy table is set to INDEF, which signals compband and 
	# compspec to compute the wavelength set from the expression.

	if (Memc[table] == EOS || streq (Memc[table], "none")) {
	    nwave = SZ_WAVE
	    call malloc (wave, nwave, TY_REAL)
	    Memr[wave] = INDEFR

	} else {
	    tp = tbtopn (wavetable, READ_ONLY, NULL)
	    nwave = tbpsta (tp, TBL_NROWS)
	    call malloc (wave, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)
	}

	# Wavelength units are always angstroms

	call malloc (units, SZ_UNIT, TY_CHAR)
	call strcpy ("angstroms", Memc[units], SZ_UNIT)

	call sfree (sp)
end
