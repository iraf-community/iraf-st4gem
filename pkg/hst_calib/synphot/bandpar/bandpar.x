#* HISTORY *
#* B.Simon	06-Oct-94	original

# BANDPAR -- Calculate the photometric parameters of a bandpass

procedure bandpar ()

#--
pointer	obsmode		# Observation mode
pointer	output		# Output table name
pointer	photlist	# List of photometric parameters to calculate
real	refwave		# Ref wavelength for calculating monochromatic flux
pointer	wavtable	# Table containing wavelength array
pointer	grftable	# Instrument graph table
pointer	cmptable	# Component name table
real	hstarea		# Total telescope area

pointer	sp
real	clgetr()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (obsmode, SZ_LINE, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (photlist, SZ_FNAME, TY_CHAR)
	call salloc (wavtable, SZ_FNAME, TY_CHAR) 
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("obsmode", Memc[obsmode], SZ_LINE)
	call clgnone ("output", Memc[output], SZ_FNAME)
	call clgstr ("photlist", Memc[photlist], SZ_FNAME)
	refwave = clgetr ("refwave")

	call clgnone ("wavetab", Memc[wavtable], SZ_FNAME)
	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr( "area" )

	# Stash the telescope area

	call put_hstarea (hstarea)

	# Substitute 1.0 for blank obsmode

	if (Memc[obsmode] == EOS)
	    call strcpy ("1.0", Memc[obsmode], SZ_LINE)

	# Call the procedure that does all the work

	call comppar (Memc[obsmode], Memc[output], Memc[photlist], refwave,
		      Memc[wavtable], Memc[grftable], Memc[cmptable])

	call sfree (sp)
end
