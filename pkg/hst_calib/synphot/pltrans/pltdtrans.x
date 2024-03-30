include	"dtrans.h"

#* HISTORY *
#* B.Simon	21-Jul-94	original

# PLTDTRANS -- Plot the data in the transformation descriptor

procedure pltdtrans (gp, mktype, dtrans)

pointer	gp		# i: graphics descriptor
char	mktype[ARB]	# i: marker type
pointer	dtrans		# i: transformation data descriptor
#--
int	neff, ieff, type
pointer	xeff, yeff

int	linemarker(), getmarker()

begin
	if (dtrans == NULL)
	    return

	# Get the arrays from the structure

	neff = TRN_NEFF(dtrans)
	xeff = TRN_XEFF(dtrans)
	yeff = TRN_YEFF(dtrans)

	# If marker is a line, plot results as a connected line
	# If it is escaped with a bang, plot the following characters
	# Otherwise, plot the corresponding symbol

	if (linemarker (mktype) == YES) {
	    call setmarker (gp, mktype)
	    call gpline (gp, Memr[xeff], Memr[yeff], neff)

	} else if (mktype[1] == '!') {
	    do ieff = 0, neff-1
		call gtext (gp, Memr[xeff+ieff], Memr[yeff+ieff], 
			    mktype[2], "h=c;v=c")

	} else {
	    type = getmarker (mktype)
	    call gpmark (gp, Memr[xeff], Memr[yeff], neff, type, 1.0, 1.0)
	}

end
