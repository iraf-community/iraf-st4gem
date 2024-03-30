# CKMISMATCH -- Check for a mismatch between some header parameters and
# task parameters.  Use the header parameters if there is a conflict.

procedure ckmismatch( datafile, xmode, xform, ymode, yform, 
	                        xdmode, xdform, ydmode, ydform )

char	datafile[ARB]		# i: Input data file name (contains header)
char	xmode[ARB]		# i: User specified xmode 
char	xform[ARB]		# i: User specified xform
char	ymode[ARB]		# i: User specified ymode
char	yform[ARB]		# i: User specified yform
char	xdmode[ARB]		# i: xmode from header
char	xdform[ARB]		# i: xform from header
char	ydmode[ARB]		# i: ymode from header
char	ydform[ARB]		# i: yform from header

bool	mismatch
bool	strne()

begin
	# Put the strings in lower case
	call strlwr( xmode )
	call strlwr( xform )
	call strlwr( ymode )
	call strlwr( yform )
	call strlwr( xdmode )
	call strlwr( xdform )
	call strlwr( ydmode )
	call strlwr( ydform )

	mismatch = ( strne( xdmode, xmode) && strne( xmode, "none" ) ) ||
	           ( strne( xdform, xform) && strne( xform, "none" ) ) ||
	           ( strne( ydmode, ymode) && strne( ymode, "none" ) ) ||
	           ( strne( ydform, yform) && strne( yform, "none" ) )

	if ( mismatch ) {
	   call printf("Mismatch between header parameter and task parameter\n")
	   call printf("in filename:  %s\n")
	      call pargstr( datafile )
	   call printf("Using header parameters:\n%s\n%s\n%s\n%s\n")
	      call pargstr( xdmode )
	      call pargstr( xdform )
	      call pargstr( ydmode )
	      call pargstr( ydform )
	}

	# copy header parameters into output strings in any case
	call strcpy( xdmode, xmode, SZ_LINE )
	call strcpy( xdform, xform, SZ_FNAME )
	call strcpy( ydmode, ymode, SZ_LINE )
	call strcpy( ydform, yform, SZ_FNAME )
end

