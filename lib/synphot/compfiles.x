include	<tbset.h>
include "libsynphot.h"

# COMPFILES -- Create a list of component file names from the component names

procedure compfiles (verbose, comptab, mxlist, mxfile, 
		     ncomp, complist, filelist)

bool	verbose				# i: diagnostic message switch
char	comptab[ARB]			# i: component lookup table name
int	mxlist				# i: max number of components in list
int	mxfile				# i: maximum length of file name
int	ncomp				# i: number of component names
char	complist[SZ_COMPID,ARB]		# i: list of instrument components
char	filelist[mxfile,ARB]		# o: list of component file names
#--
int	ic, icomp, irow, nrow
pointer	sp, tabname, compid, badnames, tp, comp, file

string	compcol   "COMPNAME"
string	filecol   "FILENAME"
string	toomany   "compfiles: too many components for list"
string	undercomp "Component names not found in lookup table"
string	header    "\nList of throughput tables:\n"
string	bodyfmt   "\t%s\n"

bool	streq()
int	tbpsta(), gstrcpy()
pointer	opnsyntab()

errchk	opnsyntab, synphoterr, syncolptr, tbegtt

begin
	# Check to see if output array is large enough 

	if (ncomp > mxlist)
	    call syninterr (toomany, ncomp)

	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (compid, SZ_COMPID, TY_CHAR)
	call salloc (badnames, SZ_LINE, TY_CHAR)

	# Initialize output array

	do icomp = 1, ncomp
	    filelist[1,icomp] = EOS

	# Open component lookup table and get column pointers

	call lastfile (comptab, Memc[tabname], SZ_FNAME)
	tp = opnsyntab (Memc[tabname])

	nrow = tbpsta (tp, TBL_NROWS)

	call syncolptr (tp, compcol, 1, comp)
	call syncolptr (tp, filecol, 2, file)

	# Find throughput filename associated with each component

	do irow = 1, nrow {
	    call tbegtt (tp, comp, irow, Memc[compid], SZ_COMPID)
	    call strfix (Memc[compid])

	    do icomp = 1, ncomp {
		if (streq (Memc[compid], complist[1,icomp])) {
		    call tbegtt (tp, file, irow, filelist[1,icomp], mxfile)
		    call strfix (filelist[1,icomp])
		}
	    }
	}

	# Print verbose messages and check for components with no file

	if (verbose)
	    call printf (header)

	ic = 0
	do icomp = 1, ncomp {
	    if (filelist[1,icomp] == EOS) {
		ic = ic + gstrcpy (complist[1,icomp], Memc[badnames+ic], 
				   SZ_LINE-(ic+1))
		Memc[badnames+ic] = ','
		ic = ic + 1
	        call printf("comp without filename = %s \n")
	        call pargstr(complist[1,icomp])

	    } else if (verbose) {
		call printf (bodyfmt)
		call pargstr (filelist[1,icomp])
	    }
	}

	# Exit with an error message if a component does not have 
	# a corresponding file

	if (ic > 0) {

	    Memc[badnames+ic] = EOS
	    call synphoterr (undercomp, Memc[badnames])

	}

	call sfree (sp)
end
