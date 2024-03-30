# ETOSF -- Use calculated standard errors to determine the number of
# significant figures to display.

procedure etosf( npar, par, parerr, digits )

int	npar		# i: Number of parameter values
real	par[ARB]	# i: Parameter values
real	parerr[ARB]	# i: Standard errors of parameters
int	digits[ARB]	# o: Digits of precision to display
#--

int	ic

begin

	do ic = 1, npar {
	   if ( par[ic] != 0. && !IS_INDEFR (parerr[ic]) )
	      digits[ic] = int( log10 ( 10000 *  abs(par[ic]/parerr[ic]) ) )
	   else if ( par[ic] == 0. )
	      digits[ic] = 4
	   else
	      digits[ic] = 10
	}
end
