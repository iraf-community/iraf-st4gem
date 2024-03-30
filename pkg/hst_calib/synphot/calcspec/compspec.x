include	"constants.h"
include "../plspec/plspec.h"
define	RADIAN	(RSUN/PC/1000.)
define	RENORM	(PI*RADIAN*FJY*RADIAN)
define  SIZE_STR 8
define  ISM 1
define  LMC 2

# COMPSPEC -- Spectrum calculator

# jul 1988 kdh @ stsci - first version
# aug 1988 kdh @ stsci - implement spectrum * passband
# feb 1989 kdh @ stsci - implement spectrum x scalar and  bb temp
# mar 1989 kdh @ stsci - fix bug in -
# apr 1989 kdh @ stsci - add hi temp column
# may 1989 kdh @ stsci - add grid list seq
# jun 1989 kdh @ stsci - ass default wavelength set so rn won't fail
#			- add pl power laws
# Nov 1989 Dave Bazell -- SPP version
# Apr 1991 DB  make sure form string is lower case when a new file is input
# Oct 1991 DB  Eliminate error on empty script string, just return if empty
#                or if script = "none" || "NONE"
# Jul 1992 DB  Add ebmvl:  LMC extinction law

procedure compspec( command, iw, grtbl, cmptbl, nwave, wave, spec, form)

char	command[ARB]	# i : script with commands and arguments
int	iw		# i : position in script
char	grtbl[ARB]	# i : graph table name
char	cmptbl[ARB]	# i : component table name
int	nwave		# i: number of wavelengths
pointer	wave		# io: pointer to wave length array
pointer	spec		# io: pointer spectrum array
char	form[ARB]	# io: form (units) of spectrum
# --

char 	form1[SIZE_STR], pending[SIZE_STR], cebmv[SIZE_STR], ctemp[SIZE_STR]
char	word[SZ_LINE], script[SZ_LINE]
real 	sp[MAXWAVE], sp1[MAXWAVE], ebmvalue, temp
real	bnu()
int 	i, nwv, iput, iget, ip, nchar
int 	strlen(), ctowrd(), ctor(), strsearch()
pointer	wv, errmsg, spoint
bool	status
bool	streq()

string	nowave		"Wavelength array empty, abort COMPSPEC"
string	toomany		"Too many wavelengths, nwave = %d, MAXWAVE = %d"
string	checksum	"Checksum error in COMPSPEC, nwave = %d, iput = %d"
string	nulcmd		""
string	notparsed	"The word '%s' could not be parsed from the spectrum string"

begin

	status = true
	# Allocate memory
	call smark( spoint )
	call salloc( wv, nwave, TY_REAL)
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

	} else if( strlen( script ) <= 0 || 
	     strsearch(script, "none") > 0 || strsearch(script,"NONE") > 0 )
	   return

	# Put string in lower case and start with no arith ops pending
	call strlwr( script )
	call strcpy(pending,"",1)

        do i=1,nwave {
	   Memr[wv+i-1] = Memr[wave+i-1]
	   sp[i] = 0.
	}
	nwv = nwave

	# Allocate memory for spec array

	call malloc( spec, nwv, TY_REAL)

	# Get next word of the command script

	repeat {
	   if( ctowrd( script, iw, word, SZ_LINE) <= 0 ) {

	      # Finish off any pending operations
	      if( strlen( pending ) > 0 )
	         call strcpy("=",word,SZ_FNAME)

	      # or copy to output array and return
	      else {

	         # Clip away any extra wavelengths
	         iput = 1
	         do iget = 1, nwv
	            if( Memr[wave+iput-1]  ==  Memr[wv+iget-1] ) {
	               Memr[spec+iput-1] = sp[iget]
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

	   if( streq(word,"=") || streq(word,"+") || streq(word,"-") )
	      call arithop( word, pending, form1, form, nwv, Memr[wv], sp1, sp)

	   # Advanced Spectrum Functions

	   # Standard Galactic ISM Extinction
	   else if( streq(word,"ebmv") ) {
	      nchar = ctowrd( script, iw, cebmv, SIZE_STR)
	      ip = 1
	      nchar = ctor( cebmv, ip, ebmvalue)
	      call extspec( nwv, Memr[wv], sp, form, ISM, ebmvalue, sp )

	   # Extinction
	   } else if( streq(word,"ebmvl") ) {
	      nchar = ctowrd( script, iw, cebmv, SIZE_STR)
	      ip = 1
	      nchar = ctor( cebmv, ip, ebmvalue)
	      call extspec( nwv, Memr[wv], sp, form, LMC, ebmvalue, sp )

	   # Multiply by constant
	   } else if( streq(word,"x") )
	      call constmul( script, form, nwv, Memr[wv], iw, sp)

	   # Multiply by passband
           else if( streq(word,"*") )
	      call pbandmul( script, grtbl, cmptbl, nwv, wv, form, iw, sp)

	   # Renormalize
	   else if( streq(word,"rn") )
	      call renorm( script, grtbl, cmptbl, form, nwv, wv, iw, sp)

	   # Analytic Spectrum Models

	   # Constant Spectrum
	   else if( streq(word, "abmag")   ||
	            streq(word, "stmag")   ||
	            streq(word, "vegamag") ||
	            streq(word, "obmag")   ||
	            streq(word, "mjy")     ||
	            streq(word, "jy")      || 
	            streq(word, "fnu")     ||
	            streq(word, "flam")    || 
	            streq(word, "photlam") ||
	            streq(word, "photnu")  ||
	            streq(word, "counts")  ) {

	      call constspec( script, iw, nwv, sp)
	      call strcpy( word, form, SIZE_STR)
	   }

	   # Blackbody Spectrum
	   else if( streq(word,"bb") ) {
	      nchar = ctowrd( script, iw, ctemp, SZ_FNAME)
	      ip = 1
	      nchar = ctor(ctemp, ip, temp)
	      do i=1,nwv
	         sp[i] = bnu( Memr[wv+i-1], temp )
	      # Renormalize bb
	      call bmulkr( sp, RENORM, sp, nwv)
	      call strcpy("jy", form, SIZE_STR)
	   }

	   # Powerlaw
	   else if( streq(word, "pl") )
	      call powerlaw( script, nwv, Memr[wv], form, iw, sp)

	   # Hydrogen Spectrum
	   else if( streq(word,"hi") ) {
	      call hspec( script, nwv, Memr[wv], form, iw, sp)
	      call bmulkr( sp, RENORM, sp, nwv )

	   # Spectrum Grid
	   }  else if( streq(word,"grid") )
	      call specgrid( script, nwv, Memr[wv], form, iw, sp)

	   # Load Spectrum from disk.  Try an ascii file
	   else if ( strsearch( word, ".dat" ) > 0 ) 
	      call ascspec( nwv, Memr[wv], sp, word, form, status )

	   # If that didn't work, try a table
	   else
	      call evalspec( nwv, Memr[wv], sp, word, "FLUX", form, status)

	   # Make sure form is lower case
	   call strlwr( form )

	   # If STATUS is FALSE then the command has not been parsed so
	   # signal an error 
	   if ( !status ) {
	      call sprintf( Memc[errmsg], SZ_LINE, notparsed )
	         call pargstr( word )
	      call error( 1, Memc[errmsg] )
	   }
	}

	call sfree(spoint)

end
