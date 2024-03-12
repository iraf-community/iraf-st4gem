###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#  Synopsis:	call sfsigpar(nfree, fpar)
#  Description:	A procedure to evaluate one sigma error bars for the fit
#		parameters.
#  Arguments:	int	nfree	- number of free parameters
#		real	fpar[ARB]	- array of free parameters
#  Returns:	Updated values for sigpar[] in common.
#  Notes:	Information shared in common blocks defined in "specfit.com".
#  History:	June 1989	Gerard Kriss
#		10/30/89	gak	Changed to full error matrix
#               04/24/96	gak	Used abs(errmat[i,j]) to avoid FPE's
#
###########################################################################

include "specfit.h"

procedure sfsigpar(nfree, fpar)
int	nfree
real	fpar[ARB]

int	i, j, k,l
real	chidel, chinu, chi0
real	errmat[MAXFREE,MAXFREE]
real	sqrt()



include "specfit.com"

begin

# Get initial Chi-square
	call chispec(nfree, fpar, chi0)
	chinu = chi0 / (nfitpts - nfree)
#	chidel = 1.1 * nfree * chinu  # Approximate Avni correction for nfree
				      # "interesting" params; OK within 10%.
#	chidel = 1.0		      # Bevington's choice
	chidel = chinu		      # Fudge poor data error estimates this way

# Zero out the errors
	for ( i = 1; i <= npar; i = i + 1) {
	    for ( j = 1; j <= npar; j = j + 1) {
		sigpar[i,j] = 0.
	    }
	}

# Get the elements of the error matrix
	call geterrmat(nfree, fpar, errmat)
		


# Calculate the errors
	for ( i = 1; i <= nfree; i = i + 1 ) {
		for ( j = 1; j <= nfree; j = j + 1 ) {

		if ( i == j ) {
		  if ( errmat[i,i] < 0. ) {

			for (k=1;k<=ncomp;k=k+1) {
                		for (l=1;l<=ncpar[comtyp[k]];l=l+1) {
                      			if (parptr[l,k] == iptr[i]) {
					call eprintf("\nError at Component %d Parameter %d Errmat = %g\n")
			  		 	call pargi(k)
						call pargi(l)
						call pargr(errmat[i,i])
					call eprintf("Number is invalid - should be greater than zero, check for negative curvature.\n") 
                        		}
                	    	}
			 }


			sigpar[iptr[i],iptr[j]] = sqrt( chidel * (-errmat[i,j]))
		  } else {
			sigpar[iptr[i],iptr[j]] = sqrt( chidel * errmat[i,j] )
		  }
		} else {
		  sigpar[iptr[i],iptr[j]] = sqrt( chidel * abs(errmat[i,j]) )
		}

		}
	}
end

