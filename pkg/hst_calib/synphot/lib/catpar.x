# CATPAR -- Read in a parameter spread over several parameters.  If param
# string begins with an '@' then read it in from a file.  Returns number of
# characters read as procedure value.

int procedure catpar( proot, npar, par )

char	proot[ARB]	# i: Root of the parameter name
int	npar		# i: Number of params parameter can be spread over
char	par[ARB]	# o: Complete parameter string

int	ic, nchar
int	itoc(), strlen(), getline()
char	ich[10], p[SZ_FNAME], temp[SZ_LINE], linebuff[SZ_LINE]
pointer	fd
pointer	open()

begin

	# Blank the parameter string on entry
	par[1] = EOS

	# Loop over params and concatinate to par array
	do ic = 1, npar {
	   call strcpy ( proot, p, SZ_FNAME )
	   nchar = itoc( ic, ich, 10 )
	   call strcat ( ich, p, SZ_FNAME )

	   call clgstr ( p, temp, SZ_FNAME )

	   # If parameter begins with an '@' read from file
	   if ( temp[1] == '@' ) {
	      fd = open( temp[2], READ_ONLY, TEXT_FILE )
	      while ( getline( fd, linebuff ) != EOF ) {

	         # Replace '\n' with ' ' or put a space at end of line
	         if ( linebuff[ strlen(linebuff) ] == '\n' )
	            linebuff[ strlen(linebuff) ] = ' '
	         else
	            call strcat( " ", linebuff, SZ_LINE)

	         # Concatenate lines into par
	         call strcat( linebuff, par, SZ_LINE )

	      }
	      return
	   }

	   # If last char is not a space or a tab, put a space to separate
	   # from next parameter string
	   if ( temp[strlen(temp)] != ' ' && temp[strlen(temp)] != '\t' )
	      call strcat( ' ', p, SZ_FNAME )
	   call strcat ( temp, par, SZ_LINE )
	}
	
	return( strlen( par ) )
end
