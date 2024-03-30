#* HISTORY *
#* B.Simon	16-May-94	Original

# PHOTOUTPUT -- Output the results from calcphot

procedure photoutput (form, func, mode1, mode2, spec, nband, v0, result, ispec)
		     

char	form[ARB]	# i: Output form
char	func[ARB]	# i: Ouput function
char	mode1[ARB]	# i: First obsmode expression
char	mode2[ARB]	# i: Second obsmode expression
char	spec[ARB]	# i: Spectral expression
int	nband		# i: Number of bandpasses
real	v0		# i: Value of variable zero
real	result		# i: Calulated result
int	ispec		# u: Index to result
#--
int	fval
pointer	sp, uform, ufunc

string	title1  "    VZERO            %s     Mode: %s\n"
string	title2  "    VZERO      (%s s^-1 hstarea^-1)\n"
string  title3  "    VZERO       %s(%s - %s)\n"
string  title4  "    VZERO      %s(%s) / %s(%s)\n"
string	title5  "    VZERO      %s(%s) - %s(%s)\n"

int	findfunc(), is_count(), is_magunit()

begin
	# Convert form to upper case for output

	call smark (sp)
	call salloc (uform, SZ_FNAME, TY_CHAR)
	call salloc (ufunc, SZ_FNAME, TY_CHAR)

	call strcpy (form, Memc[uform], SZ_FNAME)
	call strupr (Memc[uform])

	call strcpy (func, Memc[ufunc], SZ_FNAME)
	call strupr (Memc[ufunc])

	# Print title if this is the first spectrum

	if (ispec == 0) {
	    call printf ("Spectrum:  %s\n")
	    call pargstr (spec)

	    fval = findfunc (func)

	    if (nband == 1) {
		if (fval > 1) {
		    call printf (title1)
		    call pargstr (Memc[ufunc])
		    call pargstr (mode1)
		} else if (is_count (form) == NO) {
		    call printf (title1)
		    call pargstr (Memc[uform])
		    call pargstr (mode1)
		} else {
		    call printf (title2)
		    call pargstr (Memc[uform])
		}
	    } else {
		if (fval > 1) {
		    call printf (title3)
		    call pargstr (Memc[ufunc])
		    call pargstr (mode1)
		    call pargstr (mode2)
		} else if (is_magunit (form) == NO) {
		    call printf (title4)
		    call pargstr (Memc[uform])
		    call pargstr (mode1) 
		    call pargstr (Memc[uform])
		    call pargstr (mode2)
		} else {
		    call printf (title5)
		    call pargstr (Memc[uform])
		    call pargstr (mode1) 
		    call pargstr (Memc[uform])
		    call pargstr (mode2)
		}
	    }
	}

	# Print results

	call printf (" %8g           %8g\n") 
	call pargr (v0)
	call pargr (result)

	call flush (STDOUT)

	# Increment count of spectra

	ispec = ispec + 1
	call sfree (sp)

end
