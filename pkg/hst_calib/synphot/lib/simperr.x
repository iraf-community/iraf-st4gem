include	"../lib/amoebapar.h"
define	IVP	($1 + ($2-1)*nvert)

# SIMPERR -- Calculate the errors for amoeba parameters from the spread of
# values in the vertices of the simplex.

procedure simperr( par, npar, pvert, parerr )

real	par[ARB]		# i: Parameter values
int	npar			# i: Number of parameters
real	pvert[ARB]		# i: Array of vertex values
real	parerr[ARB]		# o: Parameter errors

#--
int	ic, jc, nvert

begin

	nvert = npar + 1

	do ic = 1, npar {

	   parerr[ic] = 0.

	   do jc = 1, nvert
	      parerr[ic] = parerr[ic] + ( pvert[IVP(jc,ic)] - par[ic] ) * 
	                                ( pvert[IVP(jc,ic)] - par[ic] )

	   parerr[ic] = sqrt( parerr[ic]/max(npar-1,1) )

	}
	
end
