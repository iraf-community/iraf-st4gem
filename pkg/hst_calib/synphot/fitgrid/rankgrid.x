#* HISTORY *
#* B. Simon	22-Sep-94	original

# RANKGRID -- Find the two smallest chi squared values 

procedure rankgrid (ngrid, chi2, better, best)

int	ngrid		# i: number of chi squared values
real	chi2[ARB]	# i: chi squared values
int	better		# o: index of second smallest value
int	best		# o: index of smallest value
#--
int	igrid
string	toosmall  "rankgrid: not enough spectra"

begin
	# Check number of spectra in grid

	if (ngrid < 2)
	    call printerr_int (toosmall, ngrid)

	# Set initial values of better and best

	if (chi2[1] <= chi2[2]) {
	    best = 1
	    better = 2
	} else {
	    best = 2
	    better = 1
	}

	# Find smallest next smallest values of chi2

	do igrid = 3, ngrid {
	    if (chi2[igrid] < chi2[best]) {
		best = igrid
	    } else if (chi2[igrid] < chi2[better]){
		better = igrid
	    }
	}

end


