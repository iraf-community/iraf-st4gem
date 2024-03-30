#* HISTORY *
#* B. Simon	29-Jul-94	original

# GETRESID -- Compute residual of the least squares fit

procedure getresid (ndata, model, data, weight, resid, chisq)

int	ndata		# i: number of data pointes
real	model[ARB]	# i: result of the expression after evaluation
real	data[ARB]	# i: observed data to be fit
real	weight[ARB]	# i: weights used in computing residuals
double	resid[ARB]	# o: residual of the fit
double	chisq		# o: sum of squares of the residual to the fit
#--
int	idat

begin
	chisq = 0.0
	do idat = 1, ndata {
	    resid[idat] = (data[idat] - model[idat]) * weight[idat]
	    chisq = chisq + resid[idat] * resid[idat]
	}

	chisq = chisq / double (ndata)

end
