include	"../lib/amoebapar.h"
include "../plspec/plspec.h"

# FITBAND -- Fit an parametrized obsmode string to photometric data using
# the Amoeba simplex method.

procedure fitband()

pointer	model		# i: Pointer to simulated reference spectrum
pointer	parstr		# i: pointer to parameter string for ref spectrum 
pointer pfile		# i: Photometry file or @filename
pointer	spfile		# i: Spectrophotometry or @filename
pointer targets		# i: Files corresponding to targets
pointer	fitfile		# i: Name of output file for final model
real	par[MAXPAR]	# i: Initial parameter values
real	dpar[MAXPAR]	# i: Initial 

#--
pointer	sp, data, fd, temppar
pointer open()
real	hstarea, delp, relp
real	clgetr()
int	npar, nchar, ic, ip, digits[MAXPAR]
int	itoc(), strreplace(), catpar(), scount(), ctowrd()
char	ich[SZ_PAR], pstr[SZ_PAR], pname[SZ_PAR,MAXPAR],tlist[SZ_LINE,MAXPHOT]
bool	strne(), streq()

# Include targetid common block to communicate with bandfunk.
include	"targets.h"

extern	errcleanup

string	datafmt	"s= %s p= %s"

begin

	# Post error handler
	call onerror( errcleanup )

	# Allocate memory
	call smark ( sp )
	call salloc( model, SZ_LINE, TY_CHAR )
	call salloc( parstr, SZ_LINE, TY_CHAR)
	call salloc( pfile, SZ_FNAME, TY_CHAR )
	call salloc( spfile, SZ_FNAME, TY_CHAR )
	call salloc( targets, SZ_FNAME, TY_CHAR )
	call salloc( data, SZ_LINE, TY_CHAR )
	call salloc( fitfile, SZ_FNAME, TY_CHAR)
	call salloc( temppar, SZ_FNAME, TY_CHAR )

	# Get cl parameters
	call clgstr( "spfile", Memc[spfile], SZ_FNAME )
	call clgstr( "pfile", Memc[pfile], SZ_FNAME )
	call clgstr( "targets", Memc[targets], SZ_FNAME )
	hstarea = clgetr ("area")

	call strlwr( Memc[targets] )
	call put_hstarea (hstarea)

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
	nchar = catpar( "bmodel", 2, Memc[model] )
	call strlwr( Memc[model] )

	# Get initial parameter values and number of parameters
	call mgetpar( "bpar", MAXPAR, pname, par, dpar )
	npar = scount( Memc[model], "#", SZ_LINE )
	call clgstr( "fitfile", Memc[fitfile], SZ_FNAME )

	# Fill offset if necessary
	delp = clgetr( "delp" )
	relp = clgetr( "relp" )
	do ic = 1, npar {
	   if ( dpar[ic] <= 0. )
	      dpar [ic] = sqrt( delp*delp + (relp*par[ic])**2 )
	}

	# Define data string: "s=Memc[spfile] p=Memc[pfile]"
	call sprintf( Memc[data], SZ_LINE, datafmt )
	   call pargstr( Memc[spfile] )
	   call pargstr( Memc[pfile] )
	
	# Replace literal parameter names with '#n' for the n'th parameter
	do ic = 1, npar {
	   nchar = itoc( ic, ich, SZ_PAR )
	   call strcpy( "#", pstr, 1 )
	   call strcat( ich, pstr, SZ_PAR)
	   call strput( "#", pname[1,ic], 1, SZ_FNAME )
	   if ( strreplace( pname[1,ic], pstr, Memc[model], SZ_LINE) <= 0 )  {
	      call printf("Parameter %s not found in model\n%s\n")
	         call pargstr( pname[2,ic] )
	         call pargstr( Memc[model] )
	      call sfree( sp )
	      return
	   }
	}

	call abandfit( Memc[model], Memc[data], npar, par, dpar, digits )

	# Call version of insertpar that allows you to specifiy number of 
	# digits 
	call insertpard( par, digits, Memc[model] )

	# Write out final model if there is a fitfile
	call strlwr( Memc[fitfile] )
	if ( strne( Memc[fitfile], "none" ) ) {
	   fd = open( Memc[fitfile], NEW_FILE, TEXT_FILE )
	   call putline( fd, Memc[model] )
	   call close (fd)
	}
	call printf("\nFinal:  %s\n")
	   call pargstr( Memc[model] )
	if (strne( Memc[fitfile], "none" ) ) {
	   call printf("Model written to %s.\n")
	   call pargstr( Memc[fitfile] )
	}

	# Put EOS at beginning of string to blank it
	Memc[parstr] = EOS

	# Put the final parameters into task parameters
	call mputpar( "bpar", npar, pname, par, dpar, digits )

	call sfree( sp )

end
