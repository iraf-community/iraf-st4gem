include <tbset.h>

# PUTSPEC -- Write results to output table

#* HISTORY *
#* B.Simon	17-Mar-93	original PUTCOUNTS
#* B.Simon	31-May-94	check for 'none' added; form added to input
#* V.Laidler    18-Jan-02	created PUTSPEC for thermback task
#* V.Laidler    18-Nov-04	changed column heading from THERMBACK to FLUX 


procedure putspec (output, form, graphtab, ocomptab, tcomptab, mode, 
		     nwave, wave, thermspec)

char	output[ARB]	# i: output table
char	form[ARB]	# i: output form
char	graphtab[ARB]	# i: graph table
char	ocomptab[ARB]	# i: optical component lookup table
char	tcomptab[ARB]	# i: thermal component lookup table
char	mode[ARB]	# i: observation mode
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thermspec[ARB]	# i: thermal background as a function of wavelength
#--
pointer	sp, tbname, tp, wavptr, specptr

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
	call tbcdef (tp, specptr, "FLUX", form, blank, TY_REAL, 1, 1)
	call tbpset (tp, TBL_MAXPAR, 10)
	call tbtcre (tp)

	# Write header parameters

	call lastfile (graphtab, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "GRFTABLE", Memc[tbname])

	call lastfile (ocomptab, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "OCMPTBL", Memc[tbname])

	call lastfile (tcomptab, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "TCMPTBL", Memc[tbname])

	call tbhadt (tp, "OBSMODE", mode)

	# Write columns

	call tbcptr (tp, wavptr, wave, 1, nwave)
	call tbcptr (tp, specptr, thermspec, 1, nwave)

	#Close table, free memory

	call tbtclo (tp)
	call sfree (sp)
end
