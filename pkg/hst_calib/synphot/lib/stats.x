define	FACSQ 0.8483037	# 1/(log10(e) * 2.5 )**2 

# STATS -- Calculate statistics.  Returns chisq, bias, and rms given two 
# input spectra.

procedure stats( tspec, tsig, bspec, ndat, form, chisq, bias, rms, nsum)

real	tspec[ARB]	# i: Upstairs spectrum
real	tsig[ARB]	# i: Sigma for Upstairs spectrum
real	bspec[ARB]	# i: Downstairs spectrum
int	ndat		# i: number of data points
char	form[ARB]	# i: Form of data
real	chisq		# o: Chi squared
real	bias		# o: Bias
real	rms		# o: Rms
int	nsum		# o: Number of points used in calculations

real	weight, chival, resid, sumw, chisqun, sumwun, biasun, rmsun
int	ic, ngood, nbad
int	strsearch()

bool	mag

begin

	chisq = 0.
	chisqun = 0.
	bias = 0.
	biasun = 0.
	rms = 0.
	rmsun = 0.
	sumw = 0.
	sumwun = 0.
	nsum = 0.
	nbad = 0.
	ngood = 0.

	mag = false
	if ( strsearch( form, "mag" ) > 0 || strsearch( form, "MAG" ) > 0 )
	   mag = true

	# Loop over data excluding INDEFRs in the data
	do ic = 1, ndat {

	   if ( !IS_INDEFR (tspec[ic]) && !IS_INDEFR (bspec[ic]) ) {

	      nsum = nsum + 1

	      # If tsig is acceptable, increment number of good points and
	      # use tsig to calculate the weight for this point
	      if ( !IS_INDEFR (tsig[ic]) && tsig[ic] > 0 ) {

	         ngood = ngood + 1

	         # Calculate chisq
	         chival = (tspec[ic] - bspec[ic]) / tsig[ic]
	         chisq = chisq + chival * chival

	         # Magnitudes?
	         if ( mag ) {
	            resid = tspec[ic] - bspec[ic]
	            weight = 1. / (tsig[ic] * tsig[ic])

	         # Not mags so convert.
	         } else if ( tspec[ic] > 0 && bspec[ic] > 0 ) {
	            resid = -2.5 * alog10( tspec[ic]/bspec[ic] )
	            weight = FACSQ * ( tspec[ic]/tsig[ic] ) **2

	         # Data point was negative so don't use it
	         } else 
	            weight = 0.

	         sumw = sumw + weight
	         bias = bias + weight * resid
	         rms = rms + weight * resid * resid

	      # tsig is bad so increment number of bad points and weight
	      # this point as 1.
	      } else {

	         nbad = nbad + 1
	         weight = 1.

	         # Calculate chisq
	         chival = tspec[ic] - bspec[ic]
	         chisqun = chisqun + chival * chival

	         # Magnitudes?
	         if ( mag )
	            resid = tspec[ic] - bspec[ic]

	         # Not mags 
	         else if ( tspec[ic] > 0 && bspec[ic] > 0 )
	            resid = -2.5 * alog10( tspec[ic]/bspec[ic] )

	         # Negative flux so don't use the data point
	         else 
	            weight = 0.

	         sumwun = sumwun + weight
	         biasun = biasun + weight * resid
	         rmsun = rmsun + weight * resid * resid	         
	      }
	   }
	}

	if ( nsum != ngood + nbad ) {
	   call printf("Error: nsum = %d, ngood = %d, nbad = %d\n")
	      call pargi(nsum)
	      call pargi(ngood)
	      call pargi(nbad)
	}

	# If no data points had valid errors then calculate chisqr from the
	# uniformly weighted data
	if ( nbad == nsum ) {
	   chisq = chisqun / max(1,nbad)
	   if ( sumwun > 0 ) {
	      bias = biasun / sumwun
	      rms = sqrt( rmsun / sumwun )
	   }

	# If some data had valid errors then calculate statistical quantities
	# from those alone.  Report number of good and bad points if some were
	# good and some were bad.
	} else {
	   chisq = chisq / max(1, ngood)
	   if ( sumw > 0 ) {
	      bias = bias / sumw
	      rms = sqrt( rms / sumw )
	   }
	}
	if ( ngood > 0 && nbad > 0 ) {
	   call printf("Statistics calculated from %d points with errors.\n")
	      call pargi( ngood )
	   call printf("Rejected %d points with no errors.\n")
	      call pargi( nbad )
	}
end
