include "../plspec/plspec.h"

# LOADTARGETS -- Load a targetid-spectrum correspondence list.

procedure loadtargets()

int	ic, ip, nchar
int	ctowrd()
char	tlist[SZ_FNAME,MAXPHOT]
pointer targets
bool	streq()

# include targets common block
include	"targets.h"

begin

	# allocate string memory
	call malloc( targets, SZ_FNAME, TY_CHAR )

	call clgstr( "targets", Memc[targets], SZ_FNAME )
	call strlwr( Memc[targets] )

	# If Memc[targets] == "none" then set targetid to "none".  Later
	# in bandchi2 we will use the actual file name from the photfile ids
	if ( streq( Memc[targets], "none" ) )
	   call strcpy( "none", targetid[1,1], SZ_FNAME )
	else {
	   # Load targetlist and parse:
	   # first word is targetid, second is specfile
	   call loadlist( Memc[targets], ntarget, tlist, MAXPHOT )
	   do ic = 1, ntarget {
	      ip = 1
	      nchar = ctowrd( tlist[1,ic], ip, targetid[1,ic], SZ_FNAME )
	      nchar = ctowrd( tlist[1,ic], ip, specfile[1,ic], SZ_FNAME )
	   }
	}

	call mfree( targets, TY_CHAR )
end
