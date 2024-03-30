include	"amoebapar.h"

# MULTIPAR -- Read several parameters of the same rootname into a parameter
# array

procedure multipar( proot, npar, pname, par, dpar )

char	proot[ARB]	# i: Root of the parameter name
int	npar		# i: Number of parameters to read
char 	pname[SZ_PAR,ARB]# o: Parameter names
real	par[ARB]	# o: Array of parameter values
real	dpar[ARB]	# o: Array of offsets

int	ic, nchar, ip
int	itoc(), strlen(), strtosub(), ctor(), nowhite()
char	ich[10], p[SZ_FNAME], temp[SZ_FNAME]

begin

	do ic = 1, npar {

	   # Make the parameter name, e.g. param1
	   nchar = itoc( ic, ich, 10 )
	   call strcpy( proot, p, SZ_FNAME )
	   call strcat( ich, p, SZ_FNAME )

	   # Get the parameter string
	   call clgstr( p, temp, SZ_FNAME )

	   # Get the parameter name
	   ip = 1
	   nchar = strtosub( temp, ip, "=", pname[1,ic], SZ_FNAME )
	   nchar = nowhite( pname[1,ic], pname[1,ic], SZ_FNAME )

	   # Convert as much as possible to real
	   ip = ip + 1
	   nchar = ctor( temp, ip, par[ic] )

	   # Try to convert more if not at end of string yet.  First try to
	   # convert from next char (ip): if its whitespace ctor will work.
	   # If its not whitespace skip a char and start to convert
	   if ( ip < strlen( temp ) ) {
	      if ( ctor( temp, ip, dpar[ic] ) > 0 )
	         return
	      ip = ip + 1
	      nchar = ctor( temp, ip+1, dpar[ic] )
	   }
	}

end
