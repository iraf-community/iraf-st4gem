include <tbset.h>

#* HISTORY *
#* B.Simon	25-May-93	original

# PUTBAND -- Write the bandpass to an output table

procedure putband (output, obsmode, grtbl, cmptbl, nrow, wave, 
		   thruput, thruerr, phot)

char	output[ARB]	# i: output table name
char	obsmode[ARB]	# i: observation mode
char	grtbl[ARB]	# i: graph table name
char	cmptbl[ARB]	# i: component lookup table name
int	nrow		# i: length of arrays
real	wave[ARB]	# i: wavelengths
real	thruput[ARB]	# i: throughputs
real	thruerr[ARB]	# i: throughput erros
real	phot[5]		# i: photometric quantities
#--
pointer	sp, tp, table, col[3]

string	wavecol   "WAVELENGTH"
string	waveunits "angstroms"
string	thrucol   "THROUGHPUT"
string	errcol    "ERROR"
string	blank     ""

pointer	tbtopn()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (table, SZ_FNAME, TY_CHAR)

	# Open the output table and define the columns

	tp = tbtopn (output, NEW_FILE, NULL)
	call tbpset (tp, TBL_MAXPAR, 10)

	call tbcdef (tp, col[1], wavecol, waveunits, blank, TY_REAL, 1, 1)
	call tbcdef (tp, col[2], thrucol, blank, blank, TY_REAL, 1, 1)
	call tbcdef (tp, col[3], errcol, blank, blank, TY_REAL, 1, 1)
	call tbtcre (tp)

	# Write the columns to the output table

	call tbcptr (tp, col[1], wave, 1, nrow)
	call tbcptr (tp, col[2], thruput, 1, nrow)
	call tbcptr (tp, col[3], thruerr, 1, nrow)

	# Write the header keywords to the output table

	call tbhadr (tp, "photflam", phot[1])
	call tbhadr (tp, "photzpt", phot[2])
	call tbhadr (tp, "photplam", phot[3])
	call tbhadr (tp, "photbw", phot[4])
	call tbhadr (tp, "aperarea", phot[5])

	call tbhadt (tp, "obsmode", obsmode)

	call lastfile (grtbl, Memc[table], SZ_FNAME)
	call tbhadt (tp, "grftable", Memc[table])

	call lastfile (cmptbl, Memc[table], SZ_FNAME)
	call tbhadt (tp, "cmptable", Memc[table])

	call tbtclo (tp)
	call sfree (sp)
end
