#* HISTORY *
#* B.Simon	28 Apr 94	rewrite, based on Dave Bazell's version

# CALCBAND -- Calculate a function of one or more passbands

procedure calcband ()

#--
pointer	command		# Expression to calculate
pointer	output		# Output table name
pointer	wavtable	# Table containing wavelength array
pointer	grftable	# Instrument graph table
pointer	cmptable	# Component name table
real	hstarea		# Total telescope area

pointer	sp
real	clgetr()

string	nilstr	""

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (command, SZ_LINE, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR) 
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("obsmode", Memc[command], SZ_LINE)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgnone ("wavetab", Memc[wavtable], SZ_FNAME)
	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr( "area" )

	# Stash the telescope area

	call put_hstarea (hstarea)

	# Substitute 1.0 for blank obsmode

	if (Memc[command] == EOS)
	    call strcpy ("1.0", Memc[command], SZ_LINE)

	# Call the procedure that does all the work

	call compexpr (Memc[command], Memc[output], nilstr, nilstr,
		       Memc[wavtable], Memc[grftable], Memc[cmptable])

	call sfree (sp)
end
