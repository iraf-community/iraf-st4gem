# RANKFIT -- Rank a given fit from fitgrid among other fits

procedure rankfit( chi2, ebmv, scale, ispec, bchi2, bebmv, bscale, bmodel,
	           rank, nkeep, mxkeep)

real	chi2		# i: Chisqr of current fit
real	ebmv		# i: Extinction of current fit
real	scale		# i: Scale factor for current fit
int	ispec		# i: Index of fit in the following arrays
real	bchi2[ARB]	# i: Chisqrs for other best fits
real	bebmv[ARB]	# i: Extinctions of other best fits
real	bscale[ARB]	# i: Scale factors of other best fits
int	bmodel[ARB]	# i: Indices of best fits
int	rank		# o: Rank of current fit
int	nkeep		# o: Number of fits being considered (initialize in
			#    calling routine
int	mxkeep		# 

int	ic, ig, ip

begin

	rank = 1
	if ( nkeep > 1 ) {
	   do ic = 1, nkeep
	      if ( chi2 > bchi2[ic] )
	         rank = rank + 1
	}

	if ( rank <= mxkeep ) {
	   if ( rank <= nkeep ) {
	      do ig = nkeep, rank, -1 {
	         ip = min ( mxkeep, ig + 1)
	         bchi2[ip] = bchi2[ig]
	         bebmv[ip] = bebmv[ig]
	         bscale[ip] = bscale[ig]
	         bmodel[ip] = bmodel[ig]
	      }
	   }
	   bchi2[rank] = chi2
	   bebmv[rank] = ebmv
	   bscale[rank] = scale
	   bmodel[rank] = ispec
	   nkeep = min ( mxkeep, nkeep + 1)
	}
end
