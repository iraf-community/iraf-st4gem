include <tbset.h>

# LOADSPHOT -- load spectrophotometry data

procedure loadsphot( nwave, file, wave, dat, sig, fwhm, form, mode, star )

int	nwave			# io: Maximum number of waves or number read
char	file[ARB]		# io: File name
real	wave[ARB]		#  o: Wavelength array
real	dat[ARB]		#  o: Spectral data
real	sig[ARB]		#  o: 1 sigma errors
real	fwhm[ARB]		#  o: resolution
char	form[ARB]		#  o: form of spectrum
char	mode[ARB]		#  o: observation mode
char	star[ARB]		#  o: star name

char	uscore
int	i, jc, nrow, badsig, ndat, nchar
int	strlen(), tbpsta(), strsearch(), nowhite()
bool	status
pointer	tp, wptr, fptr, sptr, fwptr, nulflg, errmsg
int	fnroot(), strldx()
pointer	tbtopn()

data	uscore		/ '_' /

string	invaldim	"Invalid array dimension: nwave <= 0."
string	nodatfile	"Could not read ascii file %s"
string	noname		"File name string is empty."
string	nowave		"Column WAVELENGTH not found in file %s."
string	noflux		"Column FLUX not found in file %s."

# feb 1989 kdh @ stsci - adapt from loadspec
# may 1989 kdh @ stsci - fix lup1 so [*.tab*]*.* not interpreted as a table.
# jul 1989 kdh @ stsci - use 10% error bars if none are provided
# Oct 1989 Dave Bazell - SPP version
# Oct 1990 DB - Divide input counts by hstarea (and obmag)
# Nov 1993 Bernie Simon - Fix computation of star name

begin

	call malloc( errmsg, SZ_LINE, TY_CHAR)

	# check input dimension

	if( nwave <= 0 )
	   call error(1, invaldim)

	call strlwr(file)

	if( strlen(file) <= 0)
	   call error(1, noname)

	# If file has .dat file type, try to read as an ascii file
	if ( strsearch( file, ".dat") > 0 ) {



	   call ascsphot( nwave, wave, file, dat, sig, form, mode,
	                  star, status )

	   # Signal error if file not read
	   if ( !status ) {
	      call sprintf( Memc[errmsg], SZ_LINE, nodatfile )
	      call pargstr( file )
	      call error( 1, Memc[errmsg] )
	   # Otherwise, return
	   } else
	      return
	}

	# Get star name from file name
	
	jc = fnroot (file, star, SZ_FNAME)
	jc = strldx (uscore, star)
	if (jc > 0)
	    star[jc] = EOS
	   
	# Read data from table

	nchar = nowhite( file, file, SZ_FNAME)
	tp = tbtopn(file, READ_ONLY, NULL)
	call tbcfnd(tp, "WAVELENGTH", wptr, 1)
	call tbcfnd(tp, "FLUX", fptr, 1)
	call tbcfnd(tp, "STATERROR", sptr, 1)
	call tbcfnd(tp, "FWHM", fwptr, 1)
	nrow = tbpsta(tp, TBL_NROWS)
	ndat = min(nwave, nrow)
	nwave = ndat
	call malloc(nulflg, ndat, TY_BOOL)

	if( wptr != 0)
	   call tbcgtr(tp, wptr, wave, Memb[nulflg], 1, ndat)
	else {
	   call sprintf( Memc[errmsg], SZ_LINE, nowave )
	      call pargstr( file )
	   call error( 1, Memc[errmsg] )
	}
	
	if( fptr != 0) {
	   call tbcgtr(tp, fptr, dat, Memb[nulflg], 1, ndat)
	   call tbcigt(fptr, TBL_COL_UNITS, form, SZ_FNAME)
	} else {
	   call sprintf( Memc[errmsg], SZ_LINE, noflux )
	      call pargstr( file )
	   call error( 1, Memc[errmsg] )
	}

	if ( sptr != 0)
	   call tbcgtr(tp, sptr, sig, Memb[nulflg], 1, ndat)
	else {
	   call printf("Column 'STATERROR' not found in %s.\n")
	      call pargstr( file )
	   call amovkr( INDEFR, sig, ndat )
	}

	if ( fwptr != 0)
	   call tbcgtr(tp, fwptr, fwhm, Memb[nulflg], 1, ndat)
	else {
	   call printf("Column 'FWHM' not found in %s.\n")
	      call pargstr( file )
	   call amovkr( INDEFR, fwhm, ndat)
	}

	# Count number of indefinite error bar values
	badsig = 0
	do i=1,nwave
	   if( sig[i] <= 0. || IS_INDEFR (sig[i]) )
	      badsig = badsig + 1

	# fake error bars if necessary

	if( badsig >= nwave ) {
	   call printf(" Faking error bars for %s.\n")
	      call pargstr(file)
	   if( strsearch( form, "mag") > 0 )
	      call specform( nwave, wave, dat, form, dat, "mjy", status)
	   call fakesig( nwave, dat, sig )
	   if( strsearch( form, "mag") > 0 )
	      call specform( nwave, wave, dat, "mjy", dat, form, status)
	}

	call mfree( errmsg, TY_CHAR)
	call mfree(nulflg, TY_BOOL)
	call tbtclo(tp)
end
