#* HISTORY *
#* B. Simon	04-Aug-94	original

# FLAGNEG -- Flag invalid errors as bad

procedure flagneg (nrow, weight, flag)

int	nrow		# i: length of wight and flag arrays
real	weight[ARB]	# i: errors associated with data
bool	flag[ARB]	# o: bad error flags
#--
int	irow

begin
	do irow = 1, nrow {
	    if (IS_INDEF(weight[irow])) {
		flag[irow] = true
	    } else  {
		flag[irow] = (weight[irow] <= 0.0)
	    }
	}

end
