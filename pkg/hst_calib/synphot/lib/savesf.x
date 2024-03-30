include <tbset.h>

# SAVESF -- Save the synthetic photometry in a table format readable by
# pltrans

procedure savesf( filename, npt, flux1, sig1, form1, mode1, flux2, sig2,
	          form2, mode2, starname, diff )

char	filename[ARB]	# i: Output file name
int	npt		# i: number of points to write
real	flux1[ARB]	# i: First column of fluxes
real	sig1[ARB]	# i: First column errors
char	form1[ARB]	# i: Form of first column
char	mode1[ARB]	# i: Mode of first column
real	flux2[ARB]	# i: Second column of fluxes
real	sig2[ARB]	# i: Second column errors
char	form2[ARB]	# i: Form of second column
char	mode2[ARB]	# i: Mode of second column
char	starname[SZ_FNAME,ARB] #i Starname associated with each point
char	diff[4]

int	ic, firstrow
int	tbpsta(), tbtacc()
bool	clgetb()
pointer	tp, f1, f2, s1, s2, id
pointer	tbtopn()

begin
	# If appending to existing file
	if ( clgetb( "append") && tbtacc(filename) == YES ) {

	   # Open table in append mode and get column pointers
	   tp = tbtopn( filename, READ_WRITE, NULL )
	   call tbcfnd( tp, "FLUX1", f1, 1)
	   call tbcfnd( tp, "STATERROR1", s1, 1)
	   call tbcfnd( tp, "FLUX2", f2, 1 )
	   call tbcfnd( tp, "STATERROR2", s2, 1)
	   call tbcfnd( tp, "TARGETID", id, 1)

	   # Start appending after existing data
	   firstrow = tbpsta( tp, TBL_NROWS ) + 1

	# Need to open a new table
	} else {

	   tp = tbtopn( filename, NEW_FILE, NULL )
	   call tbpset( tp, TBL_ROWLEN, npt )
	   call tbpset( tp, TBL_MAXCOLS, 4)

	   # Get pointers to columns of interest
	   call tbcdef( tp, f1, "FLUX1", form1, " ", TY_REAL, 1, 1)
	   call tbcdef( tp, s1, "STATERROR1", form1, " ",TY_REAL, 1, 1)
	   call tbcdef( tp, f2, "FLUX2", form2, " ", TY_REAL, 1, 1)
	   call tbcdef( tp, s2, "STATERROR2", form2, " ", TY_REAL, 1, 1)
	   call tbcdef( tp, id, "TARGETID"," ", " ", -SZ_FNAME, 1, 1)

	   # Create the table
	   call tbtcre( tp )

	   # Put header parameters
	   call tbhadt( tp, "XMODE", mode1 )
	   call tbhadt( tp, "YMODE", mode2 )
	   call tbhadt( tp, "XFORM", form1 )
	   call tbhadt( tp, "YFORM", form2 )
	   call tbhadt( tp, "DIFF", diff )

	   # Table is new so firstrow = 1
	   firstrow = 1

	}

	# Check the errors to see if INDEFs are needed
	do ic = 1, npt {
	   if ( sig1[ic] <= 0. ) 
	      sig1[ic] = INDEFR
	   if ( sig2[ic] <= 0. )
	      sig2[ic] = INDEFR
	}

	# Fill the table columns with the data
	call tbcptr( tp, f1, flux1, firstrow, firstrow + npt - 1)
	call tbcptr( tp, s1, sig1, firstrow, firstrow + npt - 1)
	call tbcptr( tp, f2, flux2, firstrow, firstrow + npt - 1)
	call tbcptr( tp, s2, sig2, firstrow, firstrow + npt - 1)
	call tbcptt( tp, id, starname, SZ_FNAME, firstrow, firstrow + npt - 1)

	# Close the table
	call tbtclo( tp )

end
