define	MAXLEN 4

# INSERTPARD -- Insert parameter values into model string for for final
# display.  This version ("d") allows the user to specify the number of
# digits printed for each parameter

procedure insertpard(par, digits, model )

real	par[ARB]		# i: Parameter values
int	digits[ARB]		# i: Number of digits to print
char	model[ARB]		# o: Model string with parameter values

#--
int	ic, npar, nchar, nval
int	scount(), itoc(), strreplace()
char	chpar[SZ_FNAME], ich[MAXLEN], inum[MAXLEN]

# Include common block for communication with amoebafit and fitfunk
include	"amoebafit.h"

begin

	npar = scount( fitmodel, "#", SZ_LINE)
	call strcpy( fitmodel, model, SZ_LINE )

	do ic = 1, npar {
	   # Turn par[ic] into a character string, printing number of digits
	   # specified in digits[] array
	   call sprintf( chpar, SZ_FNAME, "%0.*g")
	      call pargi( digits[ic] )
	      call pargr( par[ic] )

	   # Turn index into a character string like '#1'
	   nchar = itoc( ic, ich, MAXLEN)
	   call strcpy( "#", inum, 1)
	   call strcat( ich, inum, MAXLEN)

	   # Replace the string "#n" with par[n]
	   nval = strreplace( inum, chpar, model, SZ_LINE)
	}

end
