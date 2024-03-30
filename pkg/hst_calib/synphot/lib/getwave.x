include		<tbset.h>
include         "../plspec/plspec.h"

# GETWAVE -- Define a default wavelength array or read array from a table
# Procedure allocates memory for 'wave' and 'units' which must be freed by user

procedure getwave(wavetable, nwave, wave, units)

char	wavetable[ARB]	# i: Name of table containing wavelength array
int	nwave		# o: Number of wavelength generated or read
pointer	wave		# o: Wavelength array
pointer units		# o: Wavelength units

real	w1log, w2log, dwave
int	iwave, iw, ic
pointer	nulflg, tp, cptr, fd, buff, tmpwav
char	wavetab[SZ_FNAME]

int	strlen(), tbpsta(), ctowrd(), strsearch(), access(), getline()
int	ctor()
pointer	tbtacc(), tbtopn(), open()

# Jul 10 92 DB Look for ascii file before STSDAS table since tables can
#              now be ascii too.
# Aug 1992 DB check that a non table is accessed as a text file
# Aug 1992 DB changed nwavedef in ascii portion of code to MAXWAVE to
#             allow more points to be read in.

begin

	iw = 1
	# If wavetable is blank or contains 'none' then create default
	# wavelength set
	if ( (ctowrd(wavetable, iw, wavetab, strlen(wavetable)) <= 0) ||
	     ( strsearch(wavetable, "none") > 0 )  ||
	     ( strsearch(wavetable, "NONE") > 0 ) ) {
	   nwave = NWAVEDEF
	   call malloc( wave, nwave, TY_REAL)
	   call malloc( units, SZ_COLUNITS, TY_CHAR) 
	   
	   w1log = log10(MINWAVE) 
	   w2log = log10(real(MAXWAVE))
	   dwave = (w2log - w1log)/real(nwave - 1)

	   do iwave = 0,nwave-1
	      Memr[wave + iwave] = 10**(w1log + iwave * dwave)

	   call strcpy( "ANGSTROMS", Memc[units], SZ_FNAME)

   	# Now try an ascii file
	} else if ( access(wavetable,0,TEXT_FILE) == YES ) {

	   fd = open( wavetable, READ_ONLY, TEXT_FILE )
	   call malloc( buff, SZ_LINE, TY_CHAR )
	   call malloc( tmpwav, MAXWAVE, TY_REAL)	   

	   # Read each line and convert to real if possible.  This
	   # allows us to skip lines beginning with comments
	   nwave = 0
	   while ( getline(fd,Memc[buff]) != EOF && nwave < MAXWAVE ) {
	      ic = 1
	      if (ctor( Memc[buff], ic, Memr[tmpwav+nwave]) > 0 )
	         nwave = nwave + 1
	   }

	   # Copy wavelength array into a buffer of the correct size
	   call malloc( wave, nwave, TY_REAL )
	   call amovr( Memr[tmpwav], Memr[wave], nwave )

	   call malloc( units, SZ_COLUNITS, TY_CHAR) 
	   call strcpy( "ANGSTROMS", Memc[units], SZ_FNAME)

	   call mfree( buff, TY_CHAR )
	   call mfree( tmpwav, TY_REAL )
	   call close( fd )

	# Try to read from table
	} else if (tbtacc( wavetable ) == YES ) {

	   tp = tbtopn( wavetable, READ_ONLY, NULL )
	   cptr = 0
	   call tbcfnd( tp, "WAVELENGTH", cptr, 1)
	   if ( cptr == 0 ){
	      call error (1, "Column WAVELENGTH not found in wavetab" )
	   }

	   nwave = tbpsta( tp, TBL_NROWS )
	   call malloc( wave, nwave, TY_REAL)
	   call malloc( units, SZ_COLUNITS, TY_CHAR)
	   call malloc( nulflg, nwave, TY_BOOL)

	   call tbcgtr( tp, cptr, Memr[wave], Memb[nulflg], 1, nwave)
	   call mfree( nulflg, TY_BOOL)

	   call tbcigt( cptr, TBL_COL_UNITS, Memc[units], SZ_FNAME)
	   call tbtclo(tp)

	} else
	  call error( 1, "Cannot open wavelength file" )
end 
