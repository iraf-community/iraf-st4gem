###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University

#  Synopsis:    procedure sf_dbread2 (database, ansrname, nfree, fpar)
#	char	database[ARB]	# Database
#	char	ansrname[ARB]	# Answers name
 
#  Description: SF_DBREAD2 -- Read a specfit database entry that describes
#	a desired net flux interval for an emission line.
 
#  Arguments:   char    database[ARB]	# Database
#	char	ansrname[ARB]	# Answers name
 
#  Returns:	int	nfree
# 		real	fpar[ARB]
 
#  Notes:	Shares data in specfit.com
 
#  History:	July 1989	Gerard Kriss
#		June 28 1994	Grimes		Error Check dtlocate(), dtegti()
 
###########################################################################
 
include <error.h>
include <pkg/dttext.h>
include "specfit.h"
 
procedure sf_dbread2 (database, ansrname)
char	database[ARB]	# Database
char	ansrname[ARB]	# Answer's name
 
int     i, j, rec, istat
pointer dt, sp, dbfile
 
int     dtlocate(), dtgeti(), fscan()
pointer dtmap1()
char 	ferror[SZ_LINE]
 
include "specfit.com"
 
errchk  dtlocate(), dtgeti(), dtgar()
 
begin
 	sp = NULL
 
 	call smark (sp)
 	call salloc (dbfile, SZ_FNAME, TY_CHAR)
	call sprintf (Memc[dbfile], SZ_FNAME, "sf%s")
		call pargstr (ansrname)
 	dt = dtmap1 (database, Memc[dbfile], READ_ONLY)
 
	iferr ( rec = dtlocate (dt, ansrname) ) {
		call sprintf(ferror,SZ_LINE,"Error in flux file sf%s: Can't find correct begin statement!")
			call pargstr(ansrname)
		call error(1,ferror)
		}



	iferr ( nlines = dtgeti(dt, rec, "intervals") ) {
		call sprintf(ferror,SZ_LINE,"Unable to locate keyword 'intervals' in file sf%s")
		call pargstr(ansrname)
		call error(1, ferror)
		}

# Scan the database file to collect the intervals and their parameters
	for ( i = 1; i <= nlines; i = i + 1) {
		istat = fscan(DT(dt))
		call gargr(lam1[i])		# Lambda start
		call gargr(lam2[i])		# Lambda finish
		call gargi(ncont)		# Component number for contin
		call gargi(npieces[i])		# Number of fit components
		call gargwrd(linename[1,i],SZ_CNAME)	# ASCII label

	# For a multi-component feature, get the component numbers
		for ( j = 1; j <= npieces[i]; j = j + 1 ) {
			istat = fscan(DT(dt))
			call gargi(cmpnum[j,i])
		}
	}
 
	call dtunmap(dt)
	call sfree(sp)
end
