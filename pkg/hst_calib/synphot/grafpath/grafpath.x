#include "libsynphot.h"

define	MAXFILE		100
define	MAXPARAM	3
define  SZ_KEYWRD       32              # Length of keyword
define  GRF_COMPID      1               # index to optical component col in graph tab
define  GRF_KEYWRD      2               # index to keyword col in graph tab
define  GRF_THMLID      5               # index to thermal component col in grap
# above defines are also in libsynphot.h under st4gem$lib/synphot
# but convention precludes including it 


# GRFPTH -- Print tables used to compute throughput for observation mode

#* HISTORY *
#* B.Simon	24-May-93	original
#* B.Simon	13-Oct-94	added call to getnaked
#* V.Laidler    10-Mar-03	Changed calling sequence to SEARCHGRAF

procedure grfpth ()

#--
pointer	obsmode		# observation mode string
pointer	grtbl		# graph table name
pointer	cmptbl		# component lookup table name

bool	verbose
int	icomp, ncomp
pointer	sp, nparam, parlist, flist, fname, mode

data	verbose  / false /

begin
	# Allocate memory for strings and temporarary arrays

	call smark (sp)
	call salloc (obsmode, SZ_FNAME, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	call salloc (nparam, MAXPARAM, TY_INT)
	call salloc (parlist, MAXPARAM*MAXFILE, TY_REAL)
	call salloc (flist, (SZ_FNAME+1)*MAXFILE, TY_CHAR)
	call salloc (mode, SZ_LINE, TY_CHAR)

	# Read task parameters

	call clgstr ("obsmode", Memc[obsmode], SZ_FNAME)
	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)

	# Produce list of table names from the observation mode

	call inisyntab
	call getnaked (Memc[obsmode], Memc[mode], SZ_LINE)
	#Calling sequence change VGL 3/10/03 (should have been 9/28/01 like the others!)
	call searchgraf (verbose, Memc[grtbl], GRF_COMPID, Memc[cmptbl], Memc[mode],
			 MAXFILE, MAXPARAM, SZ_FNAME, ncomp, Memi[nparam], 
			 Memr[parlist], Memc[flist])
	call clssyntab

	# Print the list

	fname = flist
	do icomp = 1, ncomp {
	    call printf ("%s\n")
	    call pargstr (Memc[fname])

	    fname = fname + SZ_FNAME + 1
	}

	call sfree (sp)
end
