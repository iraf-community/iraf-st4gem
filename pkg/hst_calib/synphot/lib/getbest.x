define	IVP	$1 + ($2-1)*nvert

# GETBEST -- Load the best fit parameters from amoebafit into the model string

procedure getbest( chi2, npar, pvert, par, chi2avg, ibest, model )

real	chi2[ARB]	# i: Chisq values at each of the simplex vertices
int	npar		# i: Number of parameters
real	pvert[ARB]	# i: Array of vertex values
real	par[ARB]	# o: Array of best parameter values
real	chi2avg		# o: Average value of chisq over number of vertices
int	ibest		# o: Index of best vertex in chi2 array
char	model[ARB]	# o: Model with best parameter values inserted

int	iv, ip, nvert

begin

	nvert = npar + 1

	# Find the bestfit vertex and compute averaged chisqr
	ibest = 1
	chi2avg = chi2[1]
	do iv = 2, nvert {
	   if ( chi2[iv] < chi2[ibest] )
	      ibest = iv
	   chi2avg = chi2avg + chi2[iv]
	}
	chi2avg = chi2avg / max( 1, nvert )

	do ip = 1, npar
	   par[ip] = pvert[ IVP(ibest,ip) ]
	call insertpar(par, model )

end
