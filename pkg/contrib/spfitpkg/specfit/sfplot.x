###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#  Synopsis:	procedure sfplot(specfile)
#		char	plotfile[ARB]
#  Description:	SFPLOT opens the file specified by plotfile and
#		writes out the data and the components in ASCII format
#		to be plotted later by SMONGO.
#  Arguments:	char	plotfile[ARB]	- Input file name
#  Returns:	none
#  Notes:	Shares data in "specfit.com"
#  History:	July	1989	Gerard Kriss
###########################################################################
include	"specfit.h"

procedure sfplot(plotfile)
char	plotfile[ARB]

int	i, j, fd
int	open()
real	x, val

include	"specfit.com"

begin

# Open the file
	fd = open(plotfile, APPEND, TEXT_FILE)

# Write the spectrum name at the top followed by the number of points:
	call fprintf(fd,"%s\n%d\n")
		call pargstr(specname)
		call pargi(npts)

# Write out the wavelengths, the data, the errors, the total model, and the
# contribution of each component.
# Sample range is set to the full data set for this.
	call decode_ranges("*")
	for ( i = 0; i < npts; i = i + 1 ) {
		x = Memr[lambdas+i]
		call fspec(x, Memr[fitsp+i])
		call fprintf(fd,"%10.4f %11.4g %11.4g %11.4g")
		    call pargr(x)
		    call pargr(Memr[spectrum+i])
		    call pargr(Memr[errors+i])
		    call pargr(Memr[fitsp+i])

#	Now loop through each component and evaluate it --
	    for ( j = 1; j <= ncomp; j = j + 1 ) {
		call cspec( j, x, val)
		call fprintf(fd,"%11.4g")
		    call pargr(val)
	    }
	    call fprintf(fd,"\n")
	}

# Close the file and return
	call close(fd)

end
