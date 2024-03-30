# FORMLIST -- Return a number indicating the form of the data

int procedure formlist( form )

char	form[ARB]	# form of data
bool	streq()
pointer	errmsg

string	badform		"Unknown form: %s"
string	validforms	"Valid forms:  photlam, counts, flam,\n\
              fnu, photnu, jy, mjy, abmag,\n\
              stmag,vegamag, obmag\n\n"

begin
	call strlwr( form )

	if ( streq( form, "photlam" ) )
	   return 1

	else if ( streq( form, "photpix" ) || streq( form, "counts" ) )
	   return 3

	else if ( streq( form, "flam" ) )
	   return 5

	else if ( streq( form, "fnu" ) )
	   return 7

	else if ( streq( form, "photnu" ) )
	   return 9

	else if ( streq( form, "jy" ) )
	   return 11

	else if ( streq( form, "mjy" ) )
	   return 13

	else if ( streq( form, "abmag" ) )
	   return 15

	else if ( streq( form, "stmag" ) )
	   return 17

	else if ( streq( form, "vegamag" ) )
	   return 19

	else if ( streq( form, "obmag" ) )
	   return 21

	else {
	   call eprintf( validforms )
	   call malloc( errmsg, SZ_LINE, TY_CHAR )
	   call sprintf( Memc[errmsg], SZ_LINE, badform )
	      call pargstr( form )
	   call error ( 1, Memc[errmsg] )
	}

end
