include	<tbset.h>

# CALCBAND -- Calculate a function of one or more passbands

procedure calcband

#--
pointer	command		# Expression to calculate
pointer	output		# Output table name
pointer	grftable	# Instrument graph table
pointer	cmptable	# Component name table
pointer	wavtable	# Table containing wavelength array
#--
char	ch
int	fd, nwave,  ncols, icol, iw
pointer	sp, colname, colunit, cmdline
pointer	tp, wptr, cptr, wave, spec

char	getc()
int	open(), getline()
pointer	tbtopn()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (command, SZ_LINE, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR) 
	call salloc (colname, SZ_FNAME, TY_CHAR)
	call salloc (cmdline, SZ_LINE, TY_CHAR)

	# Read task parameters

	call clgstr ("obsmode", Memc[command], SZ_LINE)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	call clgstr ("wavetab", Memc[wavtable], SZ_FNAME)

	# Read in wave array or use default
	call makewave( Memc[wavtable], nwave, wave, colunit )

	# Determine upper bound for number of columns in output table

	if (Memc[command] != '@') {
	    fd = 0
	    ncols = 2
	} else {
	    fd = open (Memc[command+1], READ_ONLY, TEXT_FILE)
	    ncols = 1
	    while (getc (fd, ch) != EOF) {
		if (ch == '\n')
		    ncols = ncols + 1
	    }
	    call seek (fd, BOF)
	}

	# Create output table

	tp = tbtopn (Memc[output], NEW_FILE, NULL)
	call tbpset (tp, TBL_ROWLEN, ncols)
	call tbpset (tp, TBL_MAXCOLS, ncols)
	call tbcdef (tp, wptr, "WAVELENGTH", Memc[colunit], " ", 
		     TY_REAL, 1, 1)
	call tbtcre (tp)

	# Compute the spectra and store in the output table

	if (ncols == 2)
	    call strcpy ("THROUGHPUT", Memc[colname], SZ_FNAME)

	if (fd == 0) {
	    iw = 1
	    call compband (Memc[command], iw, Memc[grftable], Memc[cmptable],
			   nwave, wave, spec)

	    call tbcdef (tp, cptr, Memc[colname], " ", " ", 
	                 TY_REAL, 1, 1)
	    call tbcptr (tp, cptr, Memr[spec], 1, nwave)
	    call mfree (spec, TY_REAL)

	} else {
	    icol = 1
	    while (getline (fd, Memc[cmdline]) != EOF) {
	        iw = 1
		call compband (Memc[cmdline], iw, Memc[grftable], 
			       Memc[cmptable], nwave, wave, spec)

	        if (ncols > 2) {
		    call sprintf (Memc[colname], SZ_FNAME, "THROUGHPUT%d")
		    call pargi (icol)
		}

		call tbcdef (tp, cptr, Memc[colname], " ", " ", 
		             TY_REAL, 1, 1)
		call tbcptr (tp, cptr, Memr[spec], 1, nwave)
		call mfree (spec, TY_REAL)
		icol = icol + 1
	    }
	    call close (fd)
	}

	# Close output table and release memory

	call tbcptr (tp, wptr, Memr[wave], 1, nwave)
	call tbtclo (tp)
	call mfree (wave, TY_REAL)
	call mfree (colunit, TY_CHAR)
	call sfree (sp)
end
