###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:    procedure sfgetflux (outfd, nfree, fpar)
#		int	outfd			# Output file descriptor
#		int	nfree			# number of free params
#		real	fpar[ARB]		# Array of params
#	
#  Description: SFGETFLUX -- Finds the net emission line flux contained in
#		each wavelength interval defined in the "answers" database
#
#  Arguments:	int	outfd			# Output file descriptor
#		int	nfree			# number of free params
#		real	fpar[ARB]		# Array of params
#
#  Returns:	None
#
#  Notes:	Shares data in specfit.com
# 
#  History:	July 1989	Gerard Kriss
#		August 1989	gak	Added integration of fitted spectra
#		10/27/89	gak	Added scaling by chinu to raw data errs
#		10/30/89	gak	Added covariance errors
#		10/13/92	gak	Removed scaling by chinu for sigtf2
#					Print sig2ft instead of sig2c; re-order
#		04/16/93	gak	Add bpl option for continuum subtraction
# 
###########################################################################
 

include "specfit.h"
 
procedure sfgetflux (outfd, nfree, fpar)
int	outfd
int	nfree			# number of free params
real	fpar[ARB]		# Array of params

int	index1, index2, i1, i2, i3, i4
int	i, j, k, n
real	chi0, chinu, cerrfrac, err2c, sigp
real	alph, alph1, a0, c0, c1, val, dlambda
real	a1, a2, wbreak
real	cflux, cflux1, cflux2, tflux, netflux1, netflux2, gflux
real	sig2c, sig2c1, sig2c2, sigtf2, neterr1, signet2, neterr2, siggflux
real	alog()

include	"specfit.com"

begin

# Get Chi-square for scaling errors by reduced Chi-square
	call chispec(nfree, fpar, chi0)
	chinu = chi0 / (nfitpts - nfree)

# Find a continuum interval to get scaling for continuum subtraction error.
# If one can't be found, set scale factor = 0. and use fitted errors.
	cerrfrac = 0.
	for ( i = 1; i <= nlines && npieces[i] != 0; i = i + 1 )
		{}
	i = i - 1
	if ( npieces[i] == 0 && nlines != 0 ) {
	  tflux = 0.
	  sigtf2 = 0.
	  for ( j = 0;  (Memr[lambdas+j] < lam1[i]) && (j < npts); j = j + 1) 
	  	{ }
	  for ( n = j; Memr[lambdas+n] <= lam2[i] && n < npts; n = n + 1) {
		tflux = tflux + Memr[spectrum+n] *
			(Memr[lambdas+n+1] - Memr[lambdas+n])
		sigtf2 = sigtf2 + (Memr[errors+n] *
			(Memr[lambdas+n+1] - Memr[lambdas+n]) )**2
	  }
	  cerrfrac = sqrt( sigtf2 ) / tflux
	}

# For each emission line interval calculate the continuum flux and its error,
# and then calculate the net emission by (i) integrating the raw data, or
# (ii) summing the areas of the gaussian components.

	call fprintf(outfd,
"\n\n---------------------------------LINE FLUXES------------------------------------\n\n")
	call fprintf(outfd,
"\nFeature    Wave1   Wave2  Ncomp Total      Error  Continuum    NET     Error\n\n")

	for ( i = 1; i <= nlines; i = i + 1 ) {

	# First, the continuum
	  i1 = parptr[1,ncont]
	  i2 = parptr[2,ncont]
	  if ( comtyp[ncont] == 1 ) {		# Linear continuum
	    c0 = par0[ i1 ]
	    c1 = par0[ i2 ]
	    cflux = c0 * ( lam2[i] - lam1[i] ) + (lam2[i]**2 - lam1[i]**2)*c1/2.
	    cflux = cflux / 1000.
	    sig2c = ( ( sigpar[ i1, i1 ] * (lam2[i] - lam1[i]) )**2 +
	     (sigpar[i2, i2] * (lam2[i]**2-lam1[i]**2) / 2. )**2 ) * 1.e-6 +
	     (sigpar[i1,i2] + sigpar[i2,i1]) * (lam2[i] - lam1[i]) *
		( lam2[i]**2 - lam1[i]**2 ) / 2. / 1000.
	  }
	  else if ( comtyp[ncont] == 2 ) {	# Power law
	    a0 = par0[ i1 ]
	    alph = -par0[ i2 ]
	    alph1 = 1. + alph
	    cflux = ((lam2[i]/1000.)**alph1 - (lam1[i]/1000.)**alph1) *
	    		1000. * a0 / alph1 
	    sig2c = (sigpar[ i1, i1 ] * cflux / a0)**2 +
	     (sigpar[ i2, i2 ] * ( cflux / alph1 + a0 * 1000. *
	     ((lam1[i]/1000.)**alph1 * alog(lam1[i]/1000.) -
	     (lam2[i]/1000.)**alph1 * alog(lam2[i]/1000.)) / alph1 ) )**2 +
	     (sigpar[i1,i2] + sigpar[i2,i1])  * (cflux / a0) *
	     (( cflux / alph1 + a0 * 1000. *
	     ((lam1[i]/1000.)**alph1 * alog(lam1[i]/1000.) -
	     (lam2[i]/1000.)**alph1 * alog(lam2[i]/1000.)) / alph1 ) )
	  } else if ( comtyp[ncont] == 16 ) {	# Broken Power law
	    i3 = parptr[3,ncont]
	    i4 = parptr[4,ncont]
	    a0 = par0[ i1 ]
	    wbreak = par0[ i2 ]
	    a1 = -par0[ i3 ]
	    a2 = -par0[ i4 ]
	    if ( lam2[i] <= wbreak ) {
		call getcflux(lam1[i], lam2[i], a0, a1, wbreak, cflux)
		call getsig2c(lam1[i], lam2[i], a0, a1, wbreak, sigpar[i1,i1],
		    sigpar[i1,i3], sigpar[i3,i1], sigpar[i3,i3], cflux, sig2c)
	    } else if ( lam1[i] >= wbreak ) {
		call getcflux(lam1[i], lam2[i], a0, a2, wbreak, cflux)
		call getsig2c(lam1[i], lam2[i], a0, a2, wbreak, sigpar[i1,i1],
		    sigpar[i1,i4], sigpar[i4,i1], sigpar[i4,i4], cflux, sig2c)
	    } else {
		call getcflux(lam1[i], wbreak, a0, a1, wbreak, cflux1)
		call getsig2c(lam1[i], wbreak, a0, a1, wbreak, sigpar[i1,i1],
		    sigpar[i1,i3], sigpar[i3,i1], sigpar[i3,i3], cflux1, sig2c1)
		call getcflux(wbreak, lam2[i], a0, a2, wbreak, cflux2)
		call getsig2c(wbreak, lam2[i], a0, a2, wbreak, sigpar[i1,i1],
		    sigpar[i1,i4], sigpar[i4,i1], sigpar[i4,i4], cflux2, sig2c2)
		cflux = cflux1 + cflux2
		sig2c = sig2c1 + sig2c2
	    }
	  }

	# Next, integrate the raw data and the fitted components
	# for the total flux.
	  tflux = 0.
	  sigtf2 = 0.
	  gflux = 0.
	  for ( j = 0;  (Memr[lambdas+j] < lam1[i]) && (j < npts); j = j + 1) 
	  	{ }
	  for ( n = j; Memr[lambdas+n] <= lam2[i] && n < npts; n = n + 1) {
		tflux = tflux + Memr[spectrum+n] *
			(Memr[lambdas+n+1] - Memr[lambdas+n])
		sigtf2 = sigtf2 + (Memr[errors+n] *
			(Memr[lambdas+n+1] - Memr[lambdas+n]) )**2
	  # Sum and integrate fitted components
		for ( k = 1; k <= npieces[i]; k = k + 1 ) {
		    call cspec(cmpnum[k,i], Memr[lambdas+n], val)
		    gflux = gflux + val * (Memr[lambdas+n+1] - Memr[lambdas+n])
		}
	  }
	  dlambda = Memr[lambdas+n] - Memr[lambdas+j]

	# For npieces[i] > 0, this feature is a spectral line, so get net flux.
	# If npieces[i] == 0, this is a continuum interval, so get mean flux.
	  if ( npieces[i] > 0 ) {
	  	netflux1 = tflux - cflux
		if ( cerrfrac != 0. ) {
			err2c = chinu * (cerrfrac * cflux)**2
		} else {
			err2c = sig2c
		}
	  	neterr1 = sqrt( err2c + sigtf2 )
	  } else {
		netflux1 = tflux / dlambda
		neterr1 = sqrt( sigtf2 ) / dlambda
	  }

	# Sum the total flux of the fitted contributions
	# On 10/30/89 modified to add in co-variance contributions to errors
	  netflux2 = 0.
	  signet2 = 0.
	  for ( j = 1; j <= npieces[i]; j = j + 1) {
		index1 = cmpnum[j,i]
		index2 = parptr[1, index1]
		netflux2 = netflux2 + par0[ index2 ]
		signet2 = signet2 + sigpar[ index2, index2 ]**2
		for ( k = 1; k <= npieces[i]; k = k + 1 ) {
		    if ( k != j ) {
			sigp = sigpar[ index2, parptr[1, cmpnum[k,i]] ]
			#if ( sigp > 0 && debug ) {
			#  call eprintf("sigpar[%2d,%2d] = %e.\n")
			#    call pargi(index2)
			#    call pargi(parptr[1, cmpnum[k,i]] )
			#    call pargr(sigp)
			#  sigp = - sigp
			#}
			signet2 = signet2 + sigp
		    }
		}
	  }
	  if ( signet2 < 0. ) {
		call printf("SFGETFLUX: signet2 = %e.\n")
		call pargr(signet2)
		signet2 = - signet2
	  }
	  neterr2 = sqrt( signet2 )
	  if ( npieces[i] > 0 && netflux2 > 0.0) {
	  	siggflux = neterr2 * gflux / netflux2
	  } else {
		siggflux = 0.0
	  }

	# Again, for npieces[i] == 0, get fitted continuum mean flux
	  if ( npieces[i] == 0 ) {
		netflux2 = cflux / (lam2[i] - lam1[i])
		neterr2 = sqrt( sig2c ) / (lam2[i] - lam1[i])
	  }

	# Print the results
	call fprintf(outfd,
"%-10s%8.2f%8.2f %2d%10.3g%10.3g%10.3g%10.3g%10.3g\n")
		call pargstr(linename[1,i])
		call pargr(lam1[i])
		call pargr(lam2[i])
		call pargi(0)
		call pargr(tflux)
		call pargr( sqrt(sigtf2) )
		call pargr(cflux)
		call pargr(netflux1)
		call pargr(neterr1)
	call fprintf(outfd,
"%-10s                 %2d%10.3g%10.3g          %10.3g%10.3g\n\n")
		call pargstr(linename[1,i])
		call pargi(npieces[i])
		call pargr(gflux)
		call pargr(siggflux)
		call pargr(netflux2)
		call pargr(neterr2)
	}
end

procedure getcflux(lam1, lam2, a0, a1, wbreak, cflux)
real	lam1, lam2, a0, a1, wbreak, cflux

real	alph1

begin
	    alph1 = 1. + a1
	    cflux = ((lam2/wbreak)**alph1 - (lam1/wbreak)**alph1) *
	    		wbreak * a0 / alph1 
end


procedure getsig2c(lam1, lam2, a0, a1, wbreak, s11, s12, s21, s22, cflux, sig2c)
real	lam1, lam2, a0, a1, wbreak, s11, s12, s21, s22, cflux, sig2c

real	alph1

begin
	    alph1 = 1. + a1
	    sig2c = (s11 * cflux / a0)**2 +
	     (s22 * ( cflux / alph1 + a0 * wbreak *
	     ((lam1/wbreak)**alph1 * alog(lam1/wbreak) -
	     (lam2/wbreak)**alph1 * alog(lam2/wbreak)) / alph1 ) )**2 +
	     (s12 + s21)  * (cflux / a0) *
	     (( cflux / alph1 + a0 * wbreak *
	     ((lam1/wbreak)**alph1 * alog(lam1/wbreak) -
	     (lam2/wbreak)**alph1 * alog(lam2/wbreak)) / alph1 ) )
end
