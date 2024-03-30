include <tbset.h>

# LOADPHOT -- load photometry data

procedure loadphot( ndat, file, dat, sig, form, mode, star)

int	ndat			# io: Max number of data points or number read
char	file[ARB]		# io: Name of file
real	dat[ARB]		#  o: data values read
real	sig[ARB]		#  o: 1 sigma error
char	form[SZ_FNAME,ARB]	#  o: form of data
char	mode[SZ_LINE,ARB]	#  o: observing mode
char	star[SZ_FNAME,ARB]	#  o: star names 

string	invaldim	"Invalid dimension, ndat <= 0."
string	noname		"Hey, the file name string is blank!"
string	noflux		"DATUM or FLUX column not found in file %s."
string	nodatfile	"Unable to read ascii file %s in LOADPHOT"

int	nchar, nrow, sptr, fptr, mptr, tptr, uptr
int	strlen(), scan(), tbpsta(), nowhite(), strsearch()
pointer	tp, nulflg, errmsg
pointer	tbtopn()
bool	status

# jun 1989 kdh @ stsci - adapt loadphot from loadsphot.
# Oct 1989 Dave Bazell - SPP version
# Aug 1990 DB - Add ability to read ascii files

begin
	call malloc( errmsg, SZ_LINE, TY_CHAR)

	# Check input dimension

 	if( ndat <= 0 )
	   call error(1, invaldim)

	call strlwr(file)
	nchar = nowhite( file, file, SZ_FNAME)

	# If file has .dat file type, try to read as an ascii file
	if ( strsearch( file, ".dat") > 0 ) {
	   call ascphot( ndat, file, dat, sig, form, mode, star, status )

	   # Signal error if file not read
	   if ( !status ) {
	      call sprintf( Memc[errmsg], SZ_LINE, nodatfile )
	      call pargstr( file )
	      call error( 1, Memc[errmsg] )
	   # Otherwise, return
	   } else
	      return
	}

	# Get file name interactively if it is blank

	if ( strlen(file) <= 0 ) {
	   call printf("Enter table name for photometry.\n")
	   call printf("LOADPHOT> ")
	   nchar = scan()
	      call gargstr(file)
	}
	if( strlen(file) <= 0)
	   call error(1, noname)

	# Read data from table

	tp = tbtopn(file, READ_ONLY, NULL)
	# Look for DATUM column, if not found look for "FLUX"
	call tbcfnd(tp, "DATUM", fptr,1 )
	if ( fptr == 0 )
	   call tbcfnd(tp, "FLUX", fptr, 1)
	# Look for ERROR column, if not found look for "STATERROR"
	call tbcfnd(tp, "ERROR", sptr, 1 )
	if ( sptr == 0 )
	   call tbcfnd(tp, "STATERROR", sptr, 1)
	call tbcfnd(tp, "FORM", uptr, 1)
	call tbcfnd(tp, "OBSMODE",  mptr, 1)
	call tbcfnd(tp, "TARGETID", tptr, 1)
	nrow = tbpsta(tp, TBL_NROWS)
	ndat = min(ndat, nrow)
	call malloc(nulflg, ndat, TY_BOOL)

	# Data and form
	if( fptr != 0) {
	   call tbcgtr(tp, fptr, dat, Memb[nulflg], 1, ndat)
	} else {
	   call sprintf( Memc[errmsg], SZ_LINE, noflux )
	      call pargstr( file )
	   call error( 1, Memc[errmsg] )
	}

	# Statistical Error
	if ( sptr != 0)
	   call tbcgtr(tp, sptr, sig, Memb[nulflg], 1, ndat)
	else {
	   call printf("Column 'STATERROR' or 'ERROR' not found in %s.\n")
	      call pargstr( file )
	   call amovkr( INDEFR, dat, ndat )
	}

	# Form
	if ( uptr != 0)
	   call tbcgtt( tp, uptr, form, Memb[nulflg], SZ_FNAME, 1, ndat)
	else {
	   call printf("Column 'FORM' not found in %s.\n")
	      call pargstr( file )
	}

	# Instrument Mode
	if ( mptr != 0 )
	   call tbcgtt(tp, mptr, mode, Memb[nulflg], SZ_LINE, 1, ndat)
	else {
	   call printf("Column, 'OBSMODE' not found in %s.\n")
	      call pargstr( file )
	}

	# Star name
	if (tptr != 0)
	   call tbcgtt(tp, tptr, star, Memb[nulflg], SZ_FNAME, 1, ndat)
	else {
	   call printf("Column, 'TARGETID' not found in %s.\n")
	      call pargstr( file )
	}

	call mfree( errmsg, TY_CHAR ) 
	call mfree( nulflg, TY_BOOL )
	call tbtclo( tp )

end
