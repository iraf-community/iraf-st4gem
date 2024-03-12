###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure sf_dbread (database, fitname, nfree, fpar)
#		char	database[ARB]		# Database
#		char	fitname[ARB]		# Fit name
#		int	nfree			# number of free params
#		real	fpar[ARB]		# Array of free params
#
#  Description:	SF_DBREAD -- Read a specfit database entry.
#
#  Arguments:	char	database[ARB]		# Database
#		char	fitname[ARB]		# Fit name
#
#  Returns:	int	nfree
#		real	fpar[ARB]
#
#  Notes:	Shares data in specfit.com
#
#  History:	May	1989	Gerard Kriss
#		June 28 1994	Error Check dtlocate(), dtgeti()
#		9/31/94		Grimes - Changed cmatch to not always
#					throw you out of the program
#		6/95		Grimes - Added ability to put comments
#					in files
#		6/23/95		Kriss  - Changed output formats to fixed
#					 widths to line up fields
###########################################################################

include	<error.h>
include	<pkg/dttext.h>
include "specfit.h"

procedure sf_dbread (database, fitname, nfree, fpar)
char	database[ARB]		# Database
char	fitname[ARB]		# Fit name
int	nfree
real	fpar[ARB]

int	i, j, rec, ipar, np, istat,pos
pointer	dt, sp, dbfile, complst

int	dtlocate(), dtgeti(), fscan(),stridx(), strlen()
pointer	dtmap1()
char 	ferror[SZ_LINE],temp[SZ_LINE]

include "specfit.com"

errchk	 dtgeti(), dtgar(), dtlocate()

begin
	sp = NULL

	call smark (sp)
	call salloc (dbfile, SZ_FNAME, TY_CHAR)
	call sprintf (Memc[dbfile], SZ_FNAME, "sf%s")
	    call pargstr (fitname)

	dt = dtmap1 (database, Memc[dbfile], READ_ONLY)

	iferr ( rec = dtlocate (dt, fitname) ) {
		call sprintf(ferror, SZ_LINE,"Error in database file sf%s: Can't find correct begin statement!")
			call pargstr(fitname)
		call error(1, ferror)
		}


	iferr ( ncomp = dtgeti(dt, rec, "components") ) {
		call sprintf(ferror,SZ_LINE,"Unable to locate keyword 'components' in file sf%s")
			call pargstr(fitname)
		call error(1,ferror)
		}



	call salloc(complst, ncomp*SZ_CNAME, TY_CHAR)


# Scan the database file to collect the components and decode their names
	for ( i = 0; i < ncomp; i = i + 1) {
		istat = fscan(DT(dt))
		    call gargwrd(Memc[complst+SZ_CNAME*i], SZ_CNAME)
		call cmatch(Memc[complst+SZ_CNAME*i],comtyp[i+1],1)
		call sprintf(Memc[complst+SZ_CNAME*i], SZ_CNAME, "%s%d")
		    call pargstr(Memc[complst+SZ_CNAME*i])
		    call pargi(i+1)
	}

# For each component, read in the associated parameters and their constraints
	ipar = 1
	for ( i = 0; i < ncomp; i = i + 1) {
		np = dtgeti(dt, rec, Memc[complst+SZ_CNAME*i])		
		call strcpy("",comments[1,i+1],SZ_LINE)
		call strcpy("",temp,SZ_LINE)
		call dtgstr(dt,rec,Memc[complst+SZ_CNAME*i],temp,SZ_LINE)
		pos=stridx("#",temp)
		if (pos!=0) {
			for (j=pos; j <= strlen(temp)+2;j = j + 1) {
				comments[j-pos+1,i+1]=temp[j]
			}
		}

		for ( j = 1; j <= np; j = j + 1) {
			istat = fscan (DT(dt))
			call gargr(par0[ipar])
			call gargr(blim[ipar])
			call gargr(tlim[ipar])
			call gargr(step[ipar])
			call gargr(ptol[ipar])
			call gargi(ifix[ipar])
			#call strcpy("",comments[1,ipar],SZ_LINE)
			#call gargstr(comments[1,ipar],SZ_LINE)
			parptr[j,i+1] = ipar
			iptr[ipar] = ipar
			ipar = ipar + 1
		}
	}
	npar = ipar - 1

# Find the freely varying parameters and place them in fpar
	call setpar(nfree, fpar)

	call dtunmap (dt)
	call sfree (sp)
end

# Procedure SF_DBWRITE: to write out a solution to a new database

procedure sf_dbwrite(database,fitname)
char	database[ARB], fitname[ARB]

int	i, j
pointer	dbfile, sp, dt, keywrd
int	dtmap1()

include	"specfit.com"

begin
	sp = NULL

	call smark (sp)
	call salloc (dbfile, SZ_FNAME, TY_CHAR)
	call sprintf (Memc[dbfile], SZ_FNAME, "sf%s")
	    call pargstr (fitname)
	dt = dtmap1 (database, Memc[dbfile], APPEND)

# Header info
	call dtptime(dt)
	call dtput(dt, "begin\t%s\n")
	    call pargstr(fitname)
	call dtput(dt, "\ttask\tspecfit\n")

# Write out list of component names
	call dtput(dt, "components\t%d\n")
		call pargi(ncomp)
	for ( i = 1; i <= ncomp; i = i + 1) {
		call dtput(dt, "\t\t%s\n")
			call pargstr(compkeys[1,comtyp[i]])
	}

# Write out parameters for each component
	call salloc (keywrd, SZ_FNAME, TY_CHAR)
	for ( i = 1; i <= ncomp; i = i + 1) {
		call sprintf (Memc[keywrd], SZ_FNAME, "%s%d")
		    call pargstr (compkeys[1,comtyp[i]])
		    call pargi(i)
		call dtput(dt, "\t\t%s\t%d %s\n")
		    call pargstr(Memc[keywrd])
		    call pargi(ncpar[comtyp[i]])
		    call pargstr(comments[1,i])
		for ( j = 1; j <= ncpar[comtyp[i]]; j = j + 1) {
			call dtput(dt,
			    "\t\t    %13.7g %9.3g %9.3g %12.7g %8.3g %3d\n")
			    call pargr(par0[ parptr[j,i] ])
			    call pargr(blim[ parptr[j,i] ])
			    call pargr(tlim[ parptr[j,i] ])
			    call pargr(step[ parptr[j,i] ])
			    call pargr(ptol[ parptr[j,i] ])
			    call pargi(ifix[ parptr[j,i] ])
			    #call pargstr(comments[1,parptr[j,i]])
		}
	}

	call dtunmap (dt)
	call sfree (sp)
end

# Procedure CMATCH: to match an integer numerical component type to a
#			component keyword with a minimum match criterion.

procedure cmatch(compstr, type, where)
char	compstr[ARB]
int	type, where

int	len, nmatch, i, m
int	strlen(), strncmp()


include "specfit.com"

begin
	len = strlen(compstr)
	if ( len > SZ_CNAME )	
		len = SZ_CNAME
	nmatch = 0
	for ( i = 1; i <= NKEYS; i = i + 1) {
		m = strncmp(compstr, compkeys[1,i], len)
		if ( m == 0 ) {
			type = i
			nmatch = nmatch + 1
		}
	}
	if ( nmatch == 0 ) {
		call eprintf("%s is illegal.  Allowed types are\n")
			call pargstr(compstr)
		for ( i = 1; i <= NKEYS; i = i + 1) {
			call eprintf("%s\n")
				call pargstr(compkeys[1,i])
		}
		if ( where == 1 )
			call error(1, "CMATCH: Illegal component type")
		else {
			call eprintf("CMATCH: Illegal component type\n")
			type = -1
		}
	}
	if ( nmatch > 1 ) {
		if ( where == 1 ) {
			call error(1, "CMATCH: Ambiguous component type")
		} else {
			type = -1
			call eprintf("CMATCH: Ambiguous component type\n")
		}
	}
end
