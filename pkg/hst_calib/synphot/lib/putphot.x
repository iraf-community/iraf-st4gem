include <tbset.h>
define	SIZELINE (SZ_LINE - SZ_REAL) # For tbcdef defining a long char string

# PUTPHOT -- Write a photometry data file

procedure putphot( ndat, file, dat, sig, form, mode, star, append)

int	ndat			#  i: Number of data points
char	file[ARB]		#  i: Name of output file
real	dat[ARB]		#  i: data values read
real	sig[ARB]		#  i: 1 sigma error
char	form[SZ_FNAME,ARB]	#  i: form of data
char	mode[SZ_LINE,ARB]	#  i: observing mode
char	star[SZ_LINE,ARB]	#  i: star names 
bool	append			#  i: append to existing file?

# April 1990 Dave Bazell - First code
#--

int	sptr, fptr, optr, tptr, uptr, firstrow
int	tbpsta(), tbtacc()
pointer	tp
pointer	tbtopn()

begin

	call printf("\nWriting photometry to table %s\n")
	   call pargstr( file )
	if ( append && tbtacc(file) == YES ) {
	   # Open tabel in append mode and get column pointers
	   tp = tbtopn( file, READ_WRITE, NULL )
	   call tbcfnd( tp, "DATUM", fptr, 1 )
	   if ( fptr == 0 )
	      call tbcfnd( tp, "FLUX", fptr, 1 )
	   call tbcfnd( tp, "ERROR", sptr, 1 )
	   call tbcfnd( tp, "FORM", uptr, 1 )
	   call tbcfnd( tp, "OBSMODE", optr, 1 )
	   call tbcfnd( tp, "TARGETID", tptr, 1 )

	   # Start appending after existing data
	   firstrow = tbpsta( tp, TBL_NROWS ) + 1

	# Need to open new table
	} else {

	   # Open Table and define columns.  Columns with char string have
	   # data tpe defined as -MaximumNumberOfChars
	   tp = tbtopn(file, NEW_FILE, NULL)
	   call tbpset( tp, TBL_ROWLEN, ndat )
	   call tbpset( tp, TBL_MAXCOLS, 5 )

	   # Get pointers to columns of interest
	   call tbcdef( tp, fptr, "DATUM", " ", " ", TY_REAL, 1, 1 )
	   call tbcdef( tp, sptr, "ERROR", " ", " ", TY_REAL, 1, 1 )
	   call tbcdef( tp, uptr, "FORM", " ", " ", -SZ_FNAME, 1, 1)
	   call tbcdef( tp, optr, "OBSMODE", " ", " ", -SIZELINE, 1, 1)
	   call tbcdef( tp, tptr, "TARGETID", " ", " ", -SIZELINE, 1, 1)

	   # Create the table
	   call tbtcre( tp )

	   # Table is new so firstrow = 1
	   firstrow = 1
	}

	# Fluxes
	call tbcptr( tp, fptr, dat, firstrow, firstrow + ndat - 1)

	# Statistical errors
	call tbcptr( tp, sptr, sig, firstrow, firstrow + ndat - 1)

	# Forms
	call tbcptt( tp, uptr, form, SZ_FNAME, firstrow, firstrow + ndat - 1 )

	# Observation modes
	call tbcptt( tp, optr, mode, SZ_LINE, firstrow, firstrow + ndat - 1 )

	# Target ids
	call tbcptt( tp, tptr, star, SZ_LINE, firstrow, firstrow + ndat -1 )

	call tbtclo( tp )

end
