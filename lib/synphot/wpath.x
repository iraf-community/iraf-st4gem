# WPATH -- Concatenate file names into path string

procedure wpath (mxfile, ncomp, filelist, mxpath, path)

int	mxfile			# i: maximum length of filename
int	ncomp			# i: number of filenames
char	filelist[mxfile,ARB]	# i: list of filenames
int	mxpath			# i: maximum length of path string
char	path[ARB]		# o: path string
#--
int	ic, icomp

string	title  "The following throughput tables were used: "

int	gstrcpy()

begin
	ic = gstrcpy (title, path, mxpath) + 1

	do icomp = 1, ncomp {
	    if (icomp > 1) {
		path[ic] = ','
		ic = ic + 1
		path[ic] = ' '
		ic = ic + 1
	    }

	    ic = ic + gstrcpy (filelist[1,icomp], path[ic], mxpath-ic)
	    if (ic >= mxpath-2)
		break
	}

	path[ic] = EOS
end
