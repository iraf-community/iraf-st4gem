# MKTHRU -- Convert an ascii file or geis file to a fits throughput table

procedure mkthru ()

#--
pointer	input		# input file name list
bool	title		# is first row of file column names?
bool	verbose		# display messages when files converted?

int	code
pointer	sp, ifile, ofile, errmsg, pp

bool	clgetb()
int	fntgfnb(), errget()
pointer	fntopnb()

begin
	# Allocate memory for strings

	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (ifile, SZ_FNAME, TY_CHAR)
	call salloc (ofile, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Read task parameters

	call clgstr ("input", Memc[input], SZ_FNAME)

	title = clgetb ("title")
	verbose = clgetb ("verbose")

	# Loop over file names

	pp = fntopnb (Memc[input], YES)
	while (fntgfnb (pp, Memc[ifile], SZ_FNAME) != EOF) {

	    # Build output file name

	    call mkfitsname (Memc[ifile], Memc[ofile], SZ_FNAME)

	    # Copy table

	    iferr {
		call tabcopy (Memc[ifile], Memc[ofile], title)
		call hdrcopy (Memc[ifile], Memc[ofile])

	    } then {
		# Write error message

		code = errget (Memc[errmsg], SZ_LINE)

		call eprintf ("Cannot copy %s:\n%s\n")
		call pargstr (Memc[ifile])
		call pargstr (Memc[errmsg])

	    } else {
		# Write diagnostic message

		if (verbose) {
		    call eprintf ("%s -> %s\n")
		    call pargstr (Memc[ifile])
		    call pargstr (Memc[ofile])
		}
	    }
	}

	call fntclsb (pp)
	call sfree (sp)
end
