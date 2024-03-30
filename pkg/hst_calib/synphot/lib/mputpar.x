include	"amoebapar.h"

# MPUTPAR -- Put several parameters of the same rootname into a parameter
# array

procedure mputpar( proot, npar, pname, par, dpar, digits )

char	proot[ARB]	# i: Root of the parameter name
int	npar		# i: Number of parameters to read
char 	pname[SZ_PAR,ARB]# i: Parameter names
real	par[ARB]	# i Array of parameter values
real	dpar[ARB]	# i: Array of offsets
int	digits[ARB]	# i: Digits of precision

int	ic, nchar
int	itoc()
char	ich[10], p[SZ_FNAME], pstr[SZ_FNAME]

string	fmt	"%s = %0.*g, %g"

begin

	do ic = 1, npar {

	   # Make the parameter name, e.g. param1
	   nchar = itoc( ic, ich, 10 )
	   call strcpy( proot, p, SZ_FNAME )
	   call strcat( ich, p, SZ_FNAME )

	   # Make the parameter string.  Skip the first letter of the 
	   # parameter name since it is a '#'.
	   call sprintf( pstr, SZ_FNAME, fmt )
	      call pargstr( pname[2,ic] )
	      call pargi( digits[ic] )
	      call pargr( par[ic] )
	      call pargr( dpar[ic] )

	   # Fill the cl parameter
	   call clpstr( p, pstr, SZ_FNAME )

	}

end
