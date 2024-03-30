#* HISTORY *
#* B. Simon	04-Aug-94	original

# CVTWEIGHT -- Convert data errors to least squares weights

procedure cvtweight (ndata, weight)

int	ndata		# i: number of weight points
real	weight[ARB]	# u: error array on input, weight array on output
#--
int	idat
real	avg

string	badweight  "Error must be positive"

begin
	# Compute the average error

	avg = 0.0
	do idat = 1, ndata
	    avg = avg + weight[idat]
	avg = avg / real(ndata)

	# Convert error into weight

	do idat = 1, ndata {
	    if (weight[idat] > 0.0) {
		weight[idat] = avg / weight[idat]
	    } else {
		call printerr_real (badweight, weight[idat])
	    }
	}

end
