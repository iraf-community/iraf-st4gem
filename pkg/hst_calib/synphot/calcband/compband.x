include	<tbset.h>

# COMPBAND -- Passband calculator

define MAXWAVE 40000
define SIZE_STR 8

procedure compband( command, iw, grtbl, cmptbl, nwave, wave, outband)

char	command[ARB]	# i : script with commands and arguments
int	iw		# i : position in command script
char	grtbl[ARB]	# i : graph table name
char	cmptbl[ARB]	# i : component table name
int	nwave		# i:  number of wavelengths
pointer	wave		# io: pointer to wave length array
pointer	outband		# io: pointer passband array

char 	pending[SIZE_STR], form[TBL_COL_UNITS], notform[SZ_FNAME]
char	word[SZ_LINE], mode[SZ_LINE], script[SZ_LINE], testword[SZ_LINE]
real 	band[MAXWAVE], band1[MAXWAVE], const, mu, fwhm
real	fwhmlam(), avglam()
int 	i, nwv, iput, iget, ip, nchar, iw_old
int	strlen(), ctowrd(), ctor(), strdic()
pointer	wv, errmsg, spoint, flterr, filt
bool	status
bool	streq()

string	nulcmd		""
string	nowave		"Wavelength array empty, abort COMPBAND"
string	toomany		"Too many wavelengths, nwave = %d, MAXWAVE = %d"
string	noscript	"No script file in COMPBAND"
string	checksum	"Checksum error in COMPBAND, nwave = %d, iput = %d"

# Include the dictionaries containing commands, single letter modes and forms.
include "../lib/dictionaries.h"

# Nov 1989 Dave Bazell - SPP version
# Mar 1990 DB - Add call to strdic to check for compspec commands
# May 1990 DB - Add call to gmodestr to correctly complete mode string.
#  Move dict strings into header file.  Fix parsing loop for johnson 'v' etc.

begin

	# Allocate memory
	call smark( spoint )
	call salloc( wv, nwave, TY_REAL)
	call salloc( filt, nwave, TY_REAL)
	call salloc( flterr, nwave, TY_REAL)
	call salloc( errmsg, SZ_LINE, TY_CHAR)

	# Copy command to script so original is not corrupted
	call strcpy( command, script, SZ_LINE)

	# Trap invalid input
	if( nwave <= 0 ) 
	   call error(1, nowave)

	else if( nwave > MAXWAVE ) {
	   call sprintf(Memc[errmsg], SZ_LINE, toomany)
	      call pargi(nwave)
	      call pargi(MAXWAVE)
	   call error(1, Memc[errmsg])

	} else if( strlen( script ) <= 0 )
	   call error(1,noscript)

	# Put string in lower case and start with no arith ops pending
	call strlwr( script )
	call strcpy(pending,"",1)

        do i=1,nwave {
	   Memr[wv+i-1] = Memr[wave+i-1]
	   band[i] = 0.
	}
	nwv = nwave

	# Allocate memory for outband array

	call malloc( outband, nwv, TY_REAL)


	# Get next word of the command script

	repeat {

	   # Get next word from script; Save old pos'n if pushback is needed
	   iw_old = iw
	   nchar = ctowrd( script, iw, word, SZ_LINE)

	   # If no word parsed or if word is a compspec command, finish up.
	   # Compspec commands are listed in the SPECDICT string defined above
	   # and include all compspec commands except '*' which is passed
	   # directly to compband.  Must match more than one character of
	   # a compspec command or match an x (multiplication).
	   if( nchar <= 0 || 
	      ( strdic(word,testword,SZ_LINE,specdict) > 0 && nchar > 1 ) ||
	      ( streq( word,"x") && nchar == 1 ) ) {

	      # Push the word just parsed back on the stack
	      iw = iw_old

	      # Finish off any pending operations
	      if( strlen( pending ) > 0 )
	         call strcpy("=",word,SZ_FNAME)

	      # or copy to output array and return
	      else {
	         # Clip away any extra wavelengths
	         iput = 1
	         do iget = 1, nwv
	            if( Memr[wave+iput-1]  ==  Memr[wv+iget-1] ) {
	               Memr[outband+iput-1] = band[iget]
	               iput = iput + 1
	            }

	         if( iput != nwave+1 ) {
	            call sprintf( Memc[errmsg], SZ_LINE, checksum )
	               call pargi(nwave)
	               call pargi(iput-1)
	            call error(1, Memc[errmsg])
	         }
	         return
	      }
	   }

	   # Arithmetic Operations

	   if( streq(word,"=") || streq(word,"*") || streq(word,"/") )
	      call arithopb( word, pending, nwv, Memr[wv], band1, band)


	   # Advanced Passband Functions

	   # Multiply the passband by a polynomial with user specified
	   # mu and fwhm
	   else if ( streq(word, "poly") ) {
	      nchar = ctor( script, iw, mu)
	      nchar = ctor( script, iw, fwhm)
	      call poly( script, nwv, Memr[wv], mu, fwhm, iw, band )

	   # Tilt the band, i.e. multiply by polynomial.
	   # The array band contains the evaluated synthetic passband
	   } else if (streq(word, "tilt") ) {
	      fwhm = fwhmlam( nwv, Memr[wv], band )
	      mu = avglam( nwv, Memr[wv], band )
	      call tilt( script, nwv, Memr[wv], mu, fwhm, iw, band)


	   # Analytic Passband Models

	   # Gaussian Passband
	   } else if( streq(word, "gauss") ) {
	      call strcpy("synthetic", mode, SZ_LINE)
	      call gauss( script, iw, nwv, Memr[wv], band)
	   }

	   # Gaussian of log wavelength
	   else if( streq(word, "lgauss") ) {
	      call strcpy("synthetic", mode, SZ_LINE)
	      call lgauss( script, iw, nwv, Memr[wv], band)
	   }

	   # Square Passband
	   else if( streq(word, "box") ) {
	      call strcpy("synthetic", mode, SZ_LINE)
	      call box( script, nwave, Memr[wv], iw, band )
	   }


	   # Try to load a table, then convert to number, 
	   # otherwise get passband
	   else {

	      # word contains a file name or a mode string

	      call strcpy("synthetic", mode, SZ_LINE)
	      # try to read an ascii file.  ascspec skips the form check
	      # if form is set to none.
	      call strcpy( "none", notform, SZ_FNAME)
	      call ascspec( nwv, Memr[wv], band, word, notform, status )
	      # try to read a table
	      if ( !status )
	         call evalspec( nwv, Memr[wv], band, word, "THROUGHPUT",
	                        form, status)

	      # No table so try to convert to number
	      if ( !status ) {

	         ip = 1
	         nchar = ctor( word, ip, const)
	         if ( nchar >= 1 )
	            call amovkr( const, band, nwv)

	         # Try to evaluate passband
	         else {
	            call strcpy(word, mode, SZ_LINE)
	            call gmodestr( script, mode, iw, specdict, banddict )
	            call evalbandx( mode, nwv, Memr[wv], grtbl, cmptbl, 
				    Memr[filt], Memr[flterr])
	            call amovr(Memr[filt], band, nwv)
	         }
	      }
	   }
	}

	call sfree(spoint)

end
