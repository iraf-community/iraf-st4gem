include	"dictionaries.h"
include "../plspec/plspec.h"
define	SZ_FORM 10

# ASCSPEC -- Read a spectrum from a two column ascii file.  Reads header 
# lines in search of a form then reads two columns of numbers corresponding
# to WAVELENGTH and FLUX.  The FLUX and ERROR are resamples on the input
# wavelength grid, if one exists.  Otherwise, the data are returned along
# with the wavelength grid of the data.

procedure ascspec( nwave, wave, flux, name, form, status )

int	nwave		# i: Number of wavelength points
real	wave[ARB]	# i: Wavelength array to resample on
real	flux[ARB]	# o: Output flux array
char	name[ARB]	# i: Name of ascii file to read
char	form[ARB]	# o: Form of flux
bool	status		# o: File read ok?
# --

int	fd, ic, npt, nchar
int	access(), getline(), open(), strlen(), strdic(), ctor(), ctowrd()
int	strsearch()
pointer	buff, word, errmsg, spwave, spflux, sp

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

	# Open the file for reading
	fd = open( name, READ_ONLY, TEXT_FILE )
	
	# Blank the form string if not equal to "none", skip the header loop
	# if form = "none"
	call strlwr( form )
	if ( strsearch(form,"none") == 0 ) {
	   form[1] = EOS

	   # Read lines until EOF is reached.  This is the header loop.
	   while ( getline( fd, Memc[buff] ) != EOF ) {

	      # If form is empty then search for form in the first few lines
	      nchar = strlen(form)
	      if ( nchar == 0 ) {

	         call strlwr( Memc[buff] )
	         ic = 1
	         while ( ctowrd( Memc[buff], ic, Memc[word], SZ_LINE) != EOS ) {

	            # Look up each word in the form dictionary
	            if ( strdic( Memc[word],Memc[word],SZ_LINE,formdict ) > 0) {
	               call strcpy( Memc[word], form, SZ_FORM)
	               break
	            }
	         }
	         # If we found the form then break from main loop, otherwise
	         # get next line
	         if ( strlen( form ) > 0 )
	            break
	      }
	   }

	   # Signal error if form not found
	   if ( strlen( form ) == 0 ) {
	      call sprintf( Memc[errmsg], SZ_LINE, noform )
	         call pargstr( name )
	      call error( 1, Memc[errmsg] )
	   }
	}	# endif for form = "none"
	
	# Read lines until an EOF is reached.  Now we expect numbers
	npt = 1
	while( getline( fd, Memc[buff] ) != EOF ) {

	   ic = 1
	   # Try to convert first word to a number
	   if ( ctor( Memc[buff], ic, Memr[spwave+npt-1] ) > 0 ) {
	
	      # If first try is successful, try second word.  If both 
	      # words were converted to numbers then increment number of
	      # data points and get next line
	      if ( ctor( Memc[buff], ic, Memr[spflux+npt-1] ) > 0) {
	         npt = npt + 1
	         next
	      }
	   }
	}

	npt = npt - 1

	# Resample the spectrum on the input wavelength set if non zero
	if ( wave[nwave] - wave[1] > 0 && nwave > 0) {
	   call linterp ( form, npt, Memr[spwave], Memr[spflux], 
	                  nwave, wave, flux )

	# Otherwise copy data to output variables
	} else {
	   call amovr( Memr[spwave], wave, npt )
	   call amovr( Memr[spflux], flux, npt )
	   nwave = npt
	}

	call close( fd )
	call sfree( sp )
	
end
