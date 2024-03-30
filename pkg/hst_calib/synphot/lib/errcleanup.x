# ERRCLEANUP -- Catch some types of errors.

procedure errcleanup( status )

int	status	# i: Error status

begin

	if ( status != OK ) {
	   call printf( "\n Exception detected.\n")
	   call printf( " Try different initial parameter values.\n")
	   call flush( STDOUT )
	}

end
