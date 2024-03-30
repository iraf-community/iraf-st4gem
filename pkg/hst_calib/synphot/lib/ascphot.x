# ASCPHOT -- Read photometry data from a five column ascii file.  Routine
# will skip any comment lines. Columns are expected in the following order:
# DATUM ERROR FORM MODE TARGETID

procedure ascphot( ndat, name, data, sigma, form, mode, targetid, status )

int	ndat			# i: Number of data points
char	name[ARB]		# i: Name of ascii file to read
real	data[ARB]		# o: Output data array
real	sigma[ARB]		# o: Error in flux
char	form[SZ_FNAME,ARB]	# o: Form of flux
char	mode[SZ_LINE,ARB]	# o: Instrument observation mode
char	targetid[SZ_FNAME,ARB]	# o: Target id
bool	status			# o: File read ok?

# Jan 1990  DB  Check for sigma = INDEFR
# --

int	fd, ic, npt, nchar
int	access(), getline(), open(), ctor(), ctowrd()
pointer	buff, sp, errmsg

string	noform	"No FORM found in ascii file %s"

begin

	# Check if file can be accessed
	if ( access( name, READ_ONLY, TEXT_FILE ) > 0 )
	   status = true
	else {
	   status = false
	   return
	}

	# Allocate Memory
	call smark( sp )
	call salloc( buff, SZ_LINE, TY_CHAR )
	call salloc( errmsg, SZ_LINE, TY_CHAR)
	
	# Open the file for reading
	fd = open( name, READ_ONLY, TEXT_FILE )

	# Read lines until an EOF is reached.
	npt = 1
	while( getline( fd, Memc[buff] ) != EOF) {

	   ic = 1
	   # Try to convert first word to a number: data value
	   if ( ctor( Memc[buff], ic, data[npt] ) > 0 ) {
	
	      # If first try is successful then not a comment line, so
	      # just read next few words
	      nchar = ctor( Memc[buff], ic, sigma[npt] )
	      nchar = ctowrd( Memc[buff], ic, form[1,npt], SZ_FNAME )
	      nchar = ctowrd( Memc[buff], ic, mode[1,npt], SZ_LINE )
	      nchar = ctowrd( Memc[buff], ic, targetid[1,npt], SZ_FNAME )

	      npt = npt + 1
	   }
	}

	ndat = npt - 1

	call close( fd )
	call sfree( sp )
	
end
