###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure sflogit(nfree, fpar, chisq)
#
#  Description:	SFLOGIT outputs info from the current fit to the log files
#
#  Arguments:	int	nfree		- Number of free params
#		real	fpar[ARB]	- Array containing the free params
#  		real	chisq		- Chisquare for the current fit
#
#  Returns:	None
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	May	1989	Gerard Kriss
#		June	1995	grimes added printout of comments from 
#					input database
#				       fixed but caused when log_files=""
#					causing a sfree of something that 
#					had never been salloc'd
#		August	1995	grimes changed printout of sample to real
#				
#
###########################################################################

include	"specfit.h"

procedure sflogit(nfree, fpar, chisq)
int	nfree
real	fpar[ARB]
real	chisq

int	i, j, k
pointer	sp, keywrd

include	"specfit.com"

begin
	# Get last value of chi-square and evaluate parameter error bars
	call sfsigpar(nfree, fpar)


	for ( k = 1; k <= nlogfd; k = k + 1)
	{
	# Report Chi-square, nfitpts, and number of free params
		call fprintf(logfd[k], "Chisquare = %g for %d points and %d freely varying parameters.\n\n")
		    call pargr(chisq)
		    call pargi(nfitpts)
		    call pargi(nfree)

	# List sample regions of data included in the fit
	    call fprintf(logfd[k], "\nRegions of data included in the fit.\n")
	    for ( i = 1; i <= nsamp; i = i + 1) {
			call fprintf(logfd[k], "Sample %2d: %12.5g to %12.5g\n")
				call pargi(i)
				call pargr(sample[1,i])
				call pargr(sample[2,i])
	    }

	# Write out list of component names
		call fprintf(logfd[k], "\ncomponents\t%d\n")
			call pargi(ncomp)
		for ( i = 1; i <= ncomp; i = i + 1) {
			call fprintf(logfd[k], "\t\t%s\n")
				call pargstr(compkeys[1,comtyp[i]])
		}

	# Write out parameters for each component
		call fprintf(logfd[k],
"\nParam Value      +/-         Low Limit    Up Limit   Step Size   Tolerance  Fix?\n")
		sp = NULL
		call smark( sp )
		call salloc (keywrd, SZ_FNAME, TY_CHAR)
		for ( i = 1; i <= ncomp; i = i + 1) {
			call sprintf (Memc[keywrd], SZ_FNAME, "%s%d")
			    call pargstr (compkeys[1,comtyp[i]])
			    call pargi(i)
			call fprintf(logfd[k], "%s\t%d %s\n")
			    call pargstr(Memc[keywrd])
			    call pargi(ncpar[comtyp[i]])
			    call pargstr(comments[1,i])
			for ( j = 1; j <= ncpar[comtyp[i]]; j = j + 1) {
				call fprintf(logfd[k], "%14.8g%12.5g%12.5g%12.5g%12.5g%12.5g%5d\n")
				    call pargr(par0[ parptr[j,i] ])
				    call pargr(sigpar[parptr[j,i],parptr[j,i]])
				    call pargr(blim[ parptr[j,i] ])
				    call pargr(tlim[ parptr[j,i] ])
				    call pargr(step[ parptr[j,i] ])
				    call pargr(ptol[ parptr[j,i] ])
				    call pargi(ifix[ parptr[j,i] ])
			}
		}
		call fprintf(logfd[k], "\n")

	# Integrate data and fit over regions of interest to get line fluxes
		if ( ncont != -1 ) {
			call sfgetflux(logfd[k], nfree, fpar)
		}

	}

	if (nlogfd != 0) {
		call sfree (sp)
	}
end




