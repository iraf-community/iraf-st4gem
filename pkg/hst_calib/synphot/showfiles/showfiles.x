# SHOWFILES -- Show filenames used in calculating synphot expression

#* HISTORY *
#* B.Simon	24-May-93	original
#* B.Simon	13-Oct-94	added call to getnaked
#* B.Simon	09-Oct-95	changed to call calcfiles

procedure showfiles ()

#--
pointer	expr		# synphot expression string
pointer	grtbl		# graph table name
pointer	cmptbl		# component lookup table name

pointer	sp, pcode

begin
	# Allocate memory for strings and temporarary arrays

	call smark (sp)
	call salloc (expr, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)
	call salloc (pcode, SZ_LINE, TY_INT)

	# Read task parameters

	call clgstr ("expr", Memc[expr], SZ_FNAME)
	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)

	# Produce list of table names from the observation mode

	call inisyntab

	call expcompile (Memc[expr], Memi[pcode], SZ_LINE)
	call calcfiles (Memi[pcode], SZ_LINE, Memc[grtbl], Memc[cmptbl], 
			"STDOUT")

	call clssyntab

	call sfree (sp)
end
