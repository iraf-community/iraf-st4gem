###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure update(nfree, fpar)
#		procedure setpar(nfree, fpar)
#
#  Description:	UPDATE updates the stored parameter values in the parameter
#			common block with the current values of the freely
#			varying parameters.  UPDATE should be called before
#			evaluating the current fit.
#		SETPAR updates the values of the freely varying parameters in
#			the fpar array.  This should be called anytime the user
#			changes the values of the parameters.
#		SETLIM imposes user-defined limits on parameters by resetting
#			values to the center of the given range.  Should be
#			called at the beginning of a major iteration loop.
#
#  Arguments:	int	nfree		- Number of free params
#		real	fpar[ARB]	- Array containing the free params
#
#  Returns:	none
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	May	1989	Gerard Kriss
#		10/19/89	gak	Changed limit imposition algorithm
#		6/24/94		grimes  Changed Error Reporting
#
###########################################################################

include	"specfit.h"

procedure update(nfree, fpar)
int	nfree
real	fpar[ARB]

int	i, j

include "specfit.com"

begin
#call printf("Entered update.\n")
	# Update values of the freely varying parameters
	for ( i = 1; i <= nfree; i = i + 1) {
		par0[ iptr[i] ] = fpar[i]
			#call printf("%2d %2d\n")
			#call pargi(i)
			#call pargi(iptr[i])
	}

	# Now update those parameters which are linked to varying ones
	# The linked component is specified in the "fix" field, and the
	# current parameter is set to the linked value using the step size
	# as a ratio.

	for ( i = 1; i <= ncomp; i = i + 1) {
		for ( j = 1; j <= ncpar[comtyp[i]]; j = j + 1) {
			if ( ifix[parptr[j,i]] > 0 ) {
			  	par0[ parptr[j,i] ] = step[parptr[j,i]] *
					par0[parptr[j,ifix[parptr[j,i]]]]
			}
		}
	}
end

procedure setpar(nfree, fpar)
int	nfree
real	fpar[ARB]

int	i

include	"specfit.com"

begin

	nfree = 0
	for ( i = 1; i <= npar; i = i + 1) {
		if ( ifix[i] == 0 ) {
			nfree = nfree + 1
			iptr[nfree] = i
			fpar[nfree] = par0[i]
			}
	}
end



procedure setlim(nfree, fpar)
int	nfree
real	fpar[ARB]

int	i,j,k

include	"specfit.com"

begin
	for ( i = 1; i <= nfree; i = i + 1) {
		# Enforce limits.  If parameter out of range,
		## reset to 2 steps from that limit.
		# reset to 1 step from limit and decrease step size by 2x.
		if ( fpar[i] <= blim[ iptr[i] ] ) {
		
			for (k=1;k<=ncomp;k=k+1) {
                		for (j=1;j<=ncpar[comtyp[k]];j=j+1) {
                      			if (parptr[j,k] == iptr[i]) {
	call eprintf("Parameter %d of Component %d exceeded lower limit.\n")
			  		 	call pargi(j)
						call pargi(k) 
                        		}
                	    	}
			 }


			# fpar[i] = blim[ iptr[i] ] + step[ iptr[i] ] * 2.
			fpar[i] = blim[ iptr[i] ] + step[ iptr[i] ] * 1.01
			# step[ iptr[i] ] = step[ iptr[i] ] / 2.
		}
		if ( fpar[i] >= tlim[ iptr[i] ] ) {
			
			for (k=1;k<=ncomp;k=k+1) {
                        	for (j=1;j<=ncpar[comtyp[k]];j=j+1) {
                       			if (parptr[j,k] == iptr[i]) {
              call eprintf("Parameter %d of Component %d exceeded upper limit.\n")
              	                  		call pargi(j)
                                		call pargi(k)
                                	}
                            	}
                          }



			# fpar[i] = tlim[ iptr[i] ] - step[ iptr[i] ] * 2.
			fpar[i] = tlim[ iptr[i] ] - step[ iptr[i] ] * 1.01
			# step[ iptr[i] ] = step[ iptr[i] ] / 2.
		}
	}
	call update(nfree, fpar)
end



procedure freezepar(nfree, fpar)
int	nfree
real	fpar[ARB]

int	i,j,k

include	"specfit.com"

begin
	for ( i = 1; i <= nfree; i = i + 1) {
		# Enforce HARD limits.  If parameter out of range,
		# reset to boundary and change status to fixed.
		if ( fpar[i] <= blim[ iptr[i] ] ) {
				

                  	for (k=1;k<=ncomp;k=k+1) {
                        	for (j=1;j<=ncpar[comtyp[k]];j=j+1) {

                            		if (parptr[j,k] == iptr[i]) {
           call eprintf("\nParameter %d of Component %d frozen at lower limit.\n")
                           		     	call pargi(j)
                                		call pargi(k)
                                	}
                            	}
                         }
 
 


          		fpar[i] = blim[ iptr[i] ]
                	ifix[ iptr[i] ] = -1



		}

		if ( fpar[i] >= tlim[ iptr[i] ] ) {


      			for (k=1;k<=ncomp;k=k+1) {
                        	for (j=1;j<=ncpar[comtyp[k]];j=j+1) {
                        		if (parptr[j,k] == iptr[i]) {
       	call eprintf("Parameter %d of Component %d frozen at upper limit.\n")
			        		  call pargi(j)
                                		  call pargi(k)
                       		     	}
                            	}
                          }

                  	fpar[i] = tlim[ iptr[i] ]
                   	ifix[ iptr[i] ] = -1
		}
	}

	call update(nfree, fpar)
	call setpar(nfree, fpar)
end












