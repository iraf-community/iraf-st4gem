include "libsynphot.h"

# RDTHRUHDR -- Read a keyword from the throughput table headers

procedure rdthruhdr (obsmode, grftable, cmptable, keyword, 
		     sepstr, output, maxch)

char	obsmode[ARB]	# i: observation mode
char	grftable[ARB]	# i: graph table name
char	cmptable[ARB]	# i: component lookup table name
char	keyword[ARB]	# i: header keyword name
char	sepstr[ARB]	# i: separator string btw keyword values
char	output[ARB]	# o: list of keyword values
int	maxch		# i: maximum length of list
#--
int	ic, icomp, ncomp
pointer	sp, tp, filelist, keyval, file

int	gstrcpy()
pointer	opnsyntab()

errchk	graffiles

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (filelist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)
	call salloc (keyval, SZ_LINE, TY_CHAR)

	# Get list of component throughput tables

	call graffiles (obsmode, grftable, cmptable, SZ_FNAME, 
			MAXLIST, ncomp, Memc[filelist])

	# Read header parameter from throughput tables

	ic = 1
	file = filelist
	output[1] = EOS

	do icomp = 1, ncomp {
	    tp = opnsyntab (Memc[file])

	    ifnoerr {
		call rdtabhdt (tp, keyword, Memc[keyval], SZ_LINE)

	    } then {
		if (ic > 1)
		    ic = ic + gstrcpy (sepstr, output[ic], maxch-ic)

		ic = ic + gstrcpy (Memc[keyval], output[ic], maxch-ic)
	    }

	    file = file + SZ_FNAME + 1
	}

	call sfree (sp)
end
