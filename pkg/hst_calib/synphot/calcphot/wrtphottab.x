include <tbset.h>

#* HISTORY *
#* B.Simon	17-May-94	Original

# WRTPHOTTAB -- Write a new row to the photometric results table

procedure wrtphottab (form, mode, spec, result, irow)
				     
char	form[ARB]	# i: Output form
char	mode[ARB]	# i: Observation mode
char	spec[ARB]	# i: Spectral expression
real	result		# i: Calculated result of synphot expression
char	func[1]		# i: Output function
char	output[1]	# i: Name of photometric table
bool	append		# i: Append results to existing table?
int	irow		# u: Next row to write in table
#--
pointer	tp		# Table descriptor
pointer	datptr		# Results column
pointer	frmptr		# Output form column
pointer obsptr		# Obervation mode column
pointer	tgtptr		# Spectrum column

pointer	sp, temp, datcol

string	blank      " "
string	notappend  "Cannot append results of different type to table"

int	tbtacc(), tbpsta(), findfunc()
pointer	tbtopn()

begin
	if (tp == NULL)
	    return

	# Substitute for vzero in expression

	call smark (sp)
	call salloc (temp, SZ_LINE, TY_CHAR)

	call fillexpr (spec, Memc[temp], SZ_LINE)

	# Write results to output table

	call tbeptr (tp, datptr, irow, result)
	call tbeptt (tp, frmptr, irow, form)
	call tbeptt (tp, obsptr, irow, mode)
	call tbeptt (tp, tgtptr, irow, Memc[temp])

	# Increment row number

	irow = irow + 1
	call sfree (sp)
	return

	# Make a new, empty photometric results table

	entry mkphottab (output, append, func, irow)

	call smark (sp)
	call salloc (datcol, SZ_COLNAME, TY_CHAR)

	if (output[1] == EOS) {
	    tp = NULL

	} else {
	    if (findfunc (func) == 1) {
		call strcpy ("COUNTRATE", Memc[datcol], SZ_COLNAME)
	    } else {
		call strcpy (func, Memc[datcol], SZ_COLNAME)
		call strupr (Memc[datcol])
	    }

	    if (append && tbtacc (output) == YES) {
		tp = tbtopn (output, READ_WRITE, NULL)

		call tbcfnd (tp, Memc[datcol], datptr, 1)
		call tbcfnd (tp, "FORM", frmptr, 1)
		call tbcfnd (tp, "OBSMODE", obsptr, 1)
		call tbcfnd (tp, "TARGETID", tgtptr, 1)

		irow = tbpsta (tp, TBL_NROWS) + 1

		if (datptr == NULL)
		    call printerr_str (notappend, output)

	    } else {
		tp = tbtopn (output, NEW_FILE, NULL)

		call tbcdef (tp, datptr, Memc[datcol], blank, blank, 
			     TY_REAL, 1, 1)
		call tbcdef (tp, frmptr, "FORM", blank, blank, -10, 1, 1)
		call tbcdef (tp, obsptr, "OBSMODE", blank, blank, -140, 1, 1)
		call tbcdef (tp, tgtptr, "TARGETID", blank, blank, -140, 1, 1)

		call tbtcre (tp)
		irow = 1
	    }
	}

	call sfree (sp)
	return

	# Close the photometic results table

	entry clsphottab ()

	if (tp != NULL)
	    call tbtclo (tp)

	return
end
