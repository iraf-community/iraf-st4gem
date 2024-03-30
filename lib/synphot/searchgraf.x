include "libsynphot.h"

# SEARCHGRAF - Search graph and component tables for list of component files

procedure searchgraf (verbose, graphtab, grf_index, comptab, mode,  mxlist, mxparam, 
		      mxfile, ncomp, nparam, paramlist, filelist)

bool	verbose				# i: diagnostic message switch
char	graphtab[ARB]			# i: graph table name
int	grf_index			# i: selects the column in the table to return
char	comptab[ARB]			# i: component lookup table name
char	mode[ARB]			# i: instrument mode string
int	mxlist				# i: max number of components in list
int	mxparam				# i: max number of params per keyword
int	mxfile				# i: maximum length of file name
int	ncomp				# o: number of instrument components
int	nparam[ARB]			# o: number of parameters per keyword
real	paramlist[mxparam,ARB]		# o: parameters of each keyword
char	filelist[mxfile,ARB]		# o: list of component file names
#--
int	nmode
pointer	sp, modenum, modelist, complist, grf

pointer	mapgraf()

errchk	breakmode, mapgraf, getpath, adjparis, compfiles

#int i
begin


	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (modenum, mxlist, TY_INT)
	call salloc (modelist, mxlist*(SZ_KEYWRD+1), TY_CHAR)
	call salloc (complist, mxlist*(SZ_COMPID+1), TY_CHAR)

	# Break mode string into list of keywords and parameters

	call breakmode (mode, mxlist, mxparam, SZ_KEYWRD, nmode, 
			nparam, paramlist, Memc[modelist])


#	for ( i = 1; i<= nmode; i=i+1 ) {
#	  call eprintf(" Keyword %d has %d modes \n")
#	  call pargi(i)
#	  call pargi(nparam[i])
#	}	
#

	# Load graph table into structure

	grf = mapgraf (graphtab)


	# Use keywords to trace path thru instrument graph

	call getpath (verbose, grf, grf_index, mxlist, nmode, ncomp, 
		      Memc[modelist], Memi[modenum], Memc[complist])

	# Adjust parameters so that they are associated with
	# components and not keywords

	call adjparlis (mxparam, nmode, ncomp, nparam, 
			Memc[modelist], Memi[modenum], paramlist)

	# Get component file associated with component names

	call compfiles (verbose, comptab, mxlist, mxfile, 
			ncomp, Memc[complist], filelist)

	call sfree (sp)
end
