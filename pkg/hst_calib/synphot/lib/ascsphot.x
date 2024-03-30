include	"dictionaries.h"
include "../plspec/plspec.h"
define	SZ_FORM 10

# ASCSPHOT -- Read spectrophotometry data from a three column ascii file.
# Reads header lines in search of a form, mode and target id then reads three
# columns of numbers corresponding to WAVELENGTH, FLUX, and ERROR

procedure ascsphot( nwave, wave, name, flux, sigma, form, mode, targetid,
	            status )

int	nwave		# i: Number of wavelength points
real	wave[ARB]	# i: Wavelength array to resample on
char	name[ARB]	# i: Name of ascii file to read
real	flux[ARB]	# o: Output flux array
real	sigma[ARB]	# o: Error in flux
char	form[ARB]	# o: Form of flux
char	mode[ARB]	# o: Instrument observation mode
char	targetid[ARB]	# o: Target id
bool	status		# o: File read ok?
# --

int	fd, ic, npt, imode, itarg, ieq, nchar
int	access(), getline(), open(), strlen(), strdic(), ctor(), ctowrd()
int	strsearch()
pointer	buff, word, errmsg, spwave, spflux, sperror, sp

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
	call salloc( word, SZ_LINE,TY_CHAR )
	call salloc( errmsg, SZ_LINE, TY_CHAR)
	call salloc( spwave, MAXWAVE, TY_REAL)
	call salloc( spflux, MAXWAVE, TY_REAL)
	call salloc( sperror, MAXWAVE, TY_REAL)
	
	# Blank the form string
	form[1] = EOS

	# Open the file for reading
	fd = open( name, READ_ONLY, TEXT_FILE )

	# Read lines until EOF is reached.  This is the header loop.
	while ( getline( fd, Memc[buff] ) != EOF ) {

	   call strlwr( Memc[buff] )

	   # If form is empty then search for form in the first few lines.
	   # form must come before 'mode=' and 'targetid='
	   if ( strlen(form) == 0 ) {

	      ic = 1
	      while ( ctowrd( Memc[buff], ic, Memc[word], SZ_LINE) != EOS ) {

	         # Look up each word in the form dictionary
	         if ( strdic( Memc[word],Memc[word],SZ_LINE,formdict ) > 0) {
	            call strcpy( Memc[word], form, SZ_FORM)
	            break
	         }
	      }

	   # Once the form has been found search for mode and targetid. 
	   # Keep reading lines of the file until mode and targetid are
	   # found or EOF is reached
	   } else {
	      imode = strsearch( Memc[buff], "mode" )
	      itarg = strsearch( Memc[buff], "targetid" )
	      if ( imode > 0 ) {
	         ieq = strsearcH( Memc[buff], "=" )
	         call strcpy( Memc[buff+ieq-1], mode, SZ_LINE )
	      }
	      if ( itarg > 0 ) {
	         ieq = strsearch( Memc[buff], "=" )
	         call strcpy( Memc[buff+ieq-1], targetid, SZ_FNAME )
	      }
	      if ( strlen( mode ) > 0 && strlen( targetid ) > 0 )
	         break
	   }
	}

	# Signal error if form not found
	if ( strlen( form ) == 0 ) {
	   call sprintf( Memc[errmsg], SZ_LINE, noform )
	      call pargstr( name )
	   call error( 1, Memc[errmsg] )
	}

	# Rewind file and start search for numbers
	call seek( fd, BOF )

	# Read lines until an EOF is reached.  Now we expect numbers
	npt = 1
	while( getline( fd, Memc[buff] ) != EOF ) {

	   ic = 1
	   # Try to convert first word to a number: wavelength
	   if ( ctor( Memc[buff], ic, Memr[spwave+npt-1] ) > 0 ) {
	
	      # If first try is successful, try second word: flux.  If
	      # successful, increment npt
	      if ( ctor( Memc[buff], ic, Memr[spflux+npt-1] ) > 0 )
	         npt = npt + 1

	      # If second try is succsessful, try third work: error.
	      # Note that npt has been incremented: must subtract 2
	      nchar = ctor( Memc[buff], ic, Memr[sperror+npt-2] )
	   }
	}

	npt = npt - 1

	# Resample the data if a wavelength array already exists
	if ( wave[nwave] - wave[1] > 0 ) {

	   # Resample the spectrum on the input wavelength set
	   call linterp ( npt, Memr[spwave], Memr[spflux], 
	                  nwave, wave, flux )
	   # Resample the errors on the input wavelength set
	   call linterp ( npt, Memr[spwave], Memr[sperror], 
	                  nwave, wave, sigma )

	# Otherwise copy data to output variables
	} else {
	   call amovr ( Memr[spwave], wave, npt )
	   call amovr ( Memr[spflux], flux, npt )
	   call amovr ( Memr[sperror], sigma, npt )
	   nwave = npt
	}

	call close( fd )
	call sfree( sp )
	
end
