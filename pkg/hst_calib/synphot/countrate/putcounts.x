include <tbset.h>

# PUTCOUNTS -- Write results to output table

#* HISTORY *
#* B.Simon	17-Mar-93	original
#* B.Simon	31-May-94	check for 'none' added; form added to input

procedure putcounts (output, form, graphtab, comptab, mode, spectrum, exptime,
		     nwave, wave, flux)

char	output[ARB]	# i: output table
char	form[ARB]	# i: output form
char	graphtab[ARB]	# i: graph table
char	comptab[ARB]	# i: component lookup table
char	mode[ARB]	# i: observation mode
char	spectrum[ARB]	# i: spectral source
real	exptime		# i: exposure time (for header)
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	flux[ARB]	# i: flux as a function of wavelength
#--
pointer	sp, tbname, tp, wavptr, fluxptr

string	blank   " "
pointer	tbtopn()

begin
	# Check to see if output table should be created

	if (output[1] == EOS)
	    return

	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tbname, SZ_FNAME, TY_CHAR)

	# Create output table
 
	tp = tbtopn (output, NEW_FILE, 0)
	call tbcdef (tp, wavptr, "WAVELENGTH", "angstroms", blank, 
		     TY_REAL, 1, 1)
	call tbcdef (tp, fluxptr, "FLUX", form, blank, TY_REAL, 1, 1)
	call tbpset (tp, TBL_MAXPAR, 10)
	call tbtcre (tp)

	# Write header parameters

	call lastfile (graphtab, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "GRFTABLE", Memc[tbname])

	call lastfile (comptab, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "CMPTABLE", Memc[tbname])

	call tbhadt (tp, "OBSMODE", mode)
	call tbhadt (tp, "SPECTRUM", spectrum)
	call tbhadr (tp, "EXPTIME", exptime)

	# Write columns

	call tbcptr (tp, wavptr, wave, 1, nwave)
	call tbcptr (tp, fluxptr, flux, 1, nwave)

	#Close table, free memory

	call tbtclo (tp)
	call sfree (sp)
end
