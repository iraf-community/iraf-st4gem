#* HISTORY *
#* B. Simon	29-Jul-94	original

# RMVNULL -- Remove nulls (indefs) from an array containing them

procedure rmvnull (lendata, ndata, flag, data)

int	lendata		# i: original length of data array 
int	ndata		# o: length of data array after nulls removed
bool	flag[ARB]	# i: null flag array
real	data[ARB]	# u: data array
#--
int	idat, jdat

begin
	# Overwrite nulls with following values

	jdat = 0
	do idat = 1, lendata {
	    if (flag[idat])
		next

	    jdat = jdat + 1
	    if (jdat < idat)
		data[jdat] = data[idat]
	}

	# Set remaining values to null

	do idat = jdat+1, lendata
	    data[idat] = INDEFR

	ndata = jdat
end
