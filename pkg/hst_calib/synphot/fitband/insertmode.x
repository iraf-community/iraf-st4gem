# INSERTMODE -- Replace the string OBSMODE in a passband model with the 
# actual mode string.

procedure insertmode( obsmode, model, passband )

char	obsmode[ARB]	# i: Instrument mode of current observation (datum)
char	model[ARB]	# io: passband model, possibly with 'obsmode' keyword
char	passband[ARB]	# o: passband model with 'obsmode' replaced

#--

int	strreplace(), strsearch()
pointer	errmsg
string	truncated "String too long:\n%s"

begin

	# Check if the keyword "obsmode" is present in the model string and
	# replace it with the mode of the current observation
	call strcpy( model, passband, SZ_LINE )
	if ( strsearch( passband, "obsmode") > 0 )
	   if( strreplace("obsmode",obsmode,passband,SZ_LINE) >= SZ_LINE ) {
	      call malloc( errmsg, SZ_LINE, TY_CHAR )
	      call sprintf( Memc[errmsg], SZ_LINE, truncated )
	         call pargstr( obsmode )
	         call pargstr( model )
	      call error( 1, Memc[errmsg] )
	   }


end
