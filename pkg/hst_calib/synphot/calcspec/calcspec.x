#* HISTORY *
#* B.Simon	28 Apr 94	rewrite, based on Dave Bazell's vversion

# CALCSPEC -- Calculate a function of one or more spectra

procedure calcspec ()

#--
pointer	command		# Expression to calculate
pointer	output		# Output table name
pointer	outform		# Form of output spectrum
pointer	vzero		# Variable list
pointer	wavtable	# Table containing wavelength array
pointer	grftable	# Instrument graph table
pointer	cmptable	# Component name table
real	hstarea		# Total telescope area

pointer	sp
real	clgetr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (command, SZ_LINE, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (outform, SZ_FNAME,  TY_CHAR) 
	call salloc (vzero, SZ_FNAME, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR) 
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("spectrum", Memc[command], SZ_LINE)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgstr ("form", Memc[outform], SZ_FNAME)
	call clgstr ("vzero", Memc[vzero], SZ_FNAME)
	call clgnone ("wavetab", Memc[wavtable], SZ_FNAME)
	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr( "area" )

	# Stash the telescope area and check the output form

	call put_hstarea (hstarea)
	call chkform (Memc[outform])

	# Call the procedure that does all the work

	call compexpr (Memc[command], Memc[output], Memc[outform], Memc[vzero],
		       Memc[wavtable], Memc[grftable], Memc[cmptable])

	call sfree (sp)
end
