include	"../plspec/plspec.h"

# GETMODE -- Get the photometry mode to use for fitband

procedure getmode( file, mode )

char	file[ARB]	# i: Name of photometry file
char	mode[ARB]	# o: Mode to use for fitband
#--

int	ic, nphot, nchoices
pointer	sp, choices[MAXPHOT]
char	modelist[SZ_LINE, MAXPHOT]
bool	strne()

begin

	call smark( sp )
	do ic = 1, MAXPHOT
	   call salloc( choices[ic], SZ_LINE, TY_CHAR )

	# Load the list of mode strings
	nphot = MAXPHOT
	call loadtext( file, "OBSMODE", modelist, nphot )

	# Save the first mode for comparison
	call strcpy( modelist[1,1], Memc[choices[1]], SZ_LINE )
	nchoices = 1
	do ic = 1, nphot {
	   if ( strne( modelist[1,ic], modelist[1,1] ) )  {
	      call strcpy( modelist[1,ic], Memc[choices[1]], SZ_LINE )
	      nchoices = nchoices + 1
	   }
	}

	if ( nchoices > 1 ) {
	   call printf( "Photometry file %s has more than one mode.\n" )
	      call pargstr( file )
	   call printf( "Choose passband to model:\n" )
	   do ic = 1, nchoices {
	      call printf( "%s\n" )
	         call pargstr( Memc[choices[ic]] )
	   }
	   call scan()
	   call gargstr( mode, SZ_LINE )

	} else
	   call strcpy( Memc[choices[1]], mode, SZ_LINE )

	
end
