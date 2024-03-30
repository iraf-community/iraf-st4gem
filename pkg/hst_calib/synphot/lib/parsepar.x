include "amoebapar.h"

# PARSEPAR -- Parse a parameter string that contains a set of parameter values
# used in fitspec.

procedure parsepar( parstr, pname, pval, oval, npar, status )

char	parstr[ARB]		# i: Parameter string
char	pname[SZ_PAR,ARB]	# o: List of parameter names
real	pval[ARB]		# o: List of parameter values
real	oval[ARB]		# o: List of initial parameter offsets
int	npar			# o: Number of parameters parsed
bool	status			# o: Status is false if error occurs in parse

int	ip, nchar, max
int	strtosub(), ctor(), nowhite(), scount()
real	delp, relp
real	clgetr()

begin

	status = true

	# Get default absolute and relative offsets for parameters
	delp = clgetr( "delp" )
	relp = clgetr( "relp" )

	# Parse each parameter name then try to convert next few chars to 
	# real.  If a ';' is encountered then check for initial offset value
	# after the ';', otherwise use default.  If a name isn't parsed or
	# a value can't be converted then return false status.

	ip = 1
	max = min( MAXPAR, scount( parstr, "=", SZ_LINE) )
	for (npar=1; npar <= max; npar=npar+1) {

	   # Parse parameter name, puting a '#' as the first char
	   # If no string parsed then break
	   if ( strtosub( parstr, ip, "=", pname[2,npar], SZ_LINE ) > 0 ) {
	      pname[1,npar] = '#'
	      nchar = nowhite( pname[1,npar], pname[1,npar], SZ_PAR )
	   } else
	      break

	   # Convert next few characters to real number
	   if ( ctor( parstr, ip, pval[npar] ) > 0 ) {

	      # If next char is a semicolon then we expect offset to follow
	      if ( parstr[ip] == ';' ) {
	         ip = ip + 1
	         if ( ctor( parstr, ip, oval[npar] ) > 0 )
	            next
	         # Didn't find a number after the semicolon so use default
	         else
	            oval[npar] = sqrt( delp*delp + (relp*pval[npar])**2 )

	      # No semicolon found so use default
	      } else
	         oval[npar] = sqrt( delp*delp + (relp*pval[npar])**2 )

	   # No numerical value could be parsed so signal an error
	   } else {
	      call printf("**Can't find value for %dth parameter %s\n")
	         call pargi( npar )
	         call pargstr( pname[1,npar] )
	      status = false
	   }

	}

	npar = npar - 1

end
