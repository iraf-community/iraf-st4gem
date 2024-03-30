#* HISTORY *
#* B. Simon	22-Sep-94	original

# GETCHISQ -- Compute the weighted chi squared for two arrays

real procedure getchisq (nwave, data, weight, flux)

int	nwave		# i: length of wave, data, and weight arrays
real	data[ARB]	# i: input spectrum data
real	weight[ARB]	# i: data point weights
real	flux[ARB]	# i: flux from grid spectrum
#--
int	iwave
real	resid, chisq

begin
	chisq = 0.0

	do iwave = 1, nwave {
	    resid = (data[iwave] - flux[iwave]) * weight[iwave]
	    chisq = chisq + resid * resid
	}

	chisq = chisq / nwave
	return (chisq)
end
