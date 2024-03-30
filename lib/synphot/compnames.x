# COMPNAMES -- Create a string containg the component file names

procedure compnames (mxlist, mxfile, mxstr, ncomp, filelist, filestr)

int	mxlist				# i: max number of components in list
int	mxfile				# i: maximum length of file name
int	mxstr				# i: maximum length of output string
int	ncomp				# i: number of component names
char	filelist[mxfile,ARB]		# i: list of component file names
char	filestr[ARB]			# o: string containing filename list
#--
int	ic, icomp

int	gstrcpy()

begin
	# Concatenate filenames into a single string

	ic = 1
	do icomp = 1, ncomp {
	    ic = ic + gstrcpy (filelist[1,icomp], filestr[ic], mxstr-(ic+1))

	    filestr[ic] = ','
	    ic = ic + 1
	    filestr[ic] = ' '
	    ic = ic + 1
	}

	# Lop off trailing comma and blank

	if (ic < 3) {
	    filestr[1] = EOS
	} else {
	    filestr[ic-2] = EOS
	}

end
