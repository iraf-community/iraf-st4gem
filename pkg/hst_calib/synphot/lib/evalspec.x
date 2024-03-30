include <tbset.h>

# EVALSPEC -- resample a spectrum on a user supplied wavelength set.
#
# E.Medeiros -- FORTRAN version 5-APR-1988
# E.Medeiros -- SPP version 15-NOV-1988
# D. Bazell -- Make into subroutine, Sep 1989
# DB add check for obmag and counts and remove area dependence Oct 90
# Apr 91 DB force form to be lower case

procedure evalspec( nwave, wave, spec, name, column, form, status)

int	nwave			# i: number of wavelength points
real	wave[ARB]		# i: wavelength array to resample on
real	spec[ARB]		# o: resampled spectrum 
char	name[ARB]	 	# i: name of table to be read and resampled
char	column[ARB]		# i: name of column to be read
char	form[ARB] 		# o: form of interpolated data
bool	status

pointer spctabptr		# spectrum data table pointer 
int	tbtopn()		# table open function value
int	tbtacc()

pointer colptr[2]		# pointer to column descriptor
pointer	spcwave			# pointer to wavelength data
pointer	spcflux			# pointer to source flux data
pointer	nulflag			# pointer to data quality flags

char	colnam[SZ_COLNAME,2]	# table column names  
char	colname[SZ_COLNAME,2]	# table column names
char	colunits[SZ_COLUNITS,2]	# table column units
char	colfmt[SZ_COLFMT,2]	# table column print format

int	nchar
int	nowhite()
int	maxspcrow		# number of rows written to source table
int	maxrow			# larger of two row counts
int	tbpsta()		# column information function value
int	colnum[2]		# table column number
int	datatype[2]		# table column datatype
int	lendata[2]		# number of elements
int	lenfmt[2]		# print format byte count
int	yn
begin

	# initialize the source table column names
	call strcpy ( "WAVELENGTH", colnam[1,1], SZ_COLNAME )
	call strcpy ( column, colnam[1,2], SZ_COLNAME )

	# Check if table exists and return if not
	nchar = nowhite( name, name, SZ_FNAME )
	yn = tbtacc(name)
	if( yn == NO ) {
	   status = false
	   return
	}
	status = true
	
	# open source table in readonly mode
	spctabptr = tbtopn ( name, READ_ONLY, NULL )

	# determine the number of rows written into the spectrum and
	# user wavelength set tables, column pointers to spectrum data, 
	# and get column information
	maxspcrow = tbpsta ( spctabptr, TBL_NROWS )
	call tbcfnd ( spctabptr, colnam, colptr, 2 )
	call tbcinf ( colptr[1], colnum[1], colname[1,1], colunits[1,1], 
		colfmt[1,1], datatype[1], lendata[1], lenfmt[1] )
	call tbcinf ( colptr[2], colnum[2], colname[1,2], colunits[1,2], 
		colfmt[1,2], datatype[2], lendata[2], lenfmt[2] )

	call strcpy( colunits[1,2], form, SZ_COLUNITS )
	call strlwr( form )

	# reserve dynamic memory for input and output data
	call malloc ( spcwave, maxspcrow, TY_REAL )
	call malloc ( spcflux, maxspcrow, TY_REAL )
	maxrow = max( nwave, maxspcrow )
	call malloc ( nulflag, maxrow, TY_BOOL )

	# get the wavelength set and input flux data from source table
	call tbcgtr ( spctabptr, colptr[1], Memr[spcwave], 
					Memb[nulflag], 1, maxspcrow )
	call tbcgtr ( spctabptr, colptr[2], Memr[spcflux], 
					Memb[nulflag], 1, maxspcrow )

	# resample the spectrum
	call linterp ( form, maxspcrow, Memr[spcwave], Memr[spcflux], 
	               nwave, wave, spec )

	call tbtclo ( spctabptr )

	call mfree ( spcwave, TY_REAL )
	call mfree ( spcflux, TY_REAL )
	call mfree ( nulflag, TY_BOOL )
end
