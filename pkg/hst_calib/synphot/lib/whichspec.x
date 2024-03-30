# WHICHSPEC -- Given a list of spectra, return an index to a specific 
# spectrum matching a given spectrum id

int procedure whichspec( specname, ic,  nspec )

char	specname[SZ_FNAME,ARB]	# i: list of spectrum names
int	ic			# i: index to specname array
int	nspec			# i: number of spectra

# include targets common block
include	"targets.h"

# Dec 1990 Dave Bazell  First Code
#--

bool	streq()
int	isp

begin

	if ( streq( targetid[1,1], "none" ) ) {
	   isp = ic
	   call strcpy( specname[1,isp], specfile[1,isp], SZ_FNAME )
	} else
	   call findsphot( specname[1,ic], targetid, ntarget, isp )

	return isp

end

