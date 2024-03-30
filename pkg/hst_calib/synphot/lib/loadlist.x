# LOADLIST -- Read in a list from an ascii file

procedure loadlist( filename, nlist, list, maxlist )

char	filename[ARB]		# i: Name of file to read
int	nlist			# o: Number of entries read from list
char	list[SZ_LINE,ARB]	# o: Array of items
int	maxlist			# i: Maximum number of items in list

int	jc, ic, item
int	getline(), open(), strsearch(), strlen()
char	linebuff[SZ_LINE]

pointer	fd, errmsg

string	noentries	"No entries in file %s"
string	toolong		"Too many items in file %s. Maximum is %d."

begin
	if ( strsearch(filename, "none") > 0 || 
	     strsearch(filename, "NONE") > 0 ) {
	   nlist = 0
	   return
	}

	call malloc(errmsg, SZ_LINE, TY_CHAR)

	fd = 0
	if ( filename[1] == '@')
	   fd = open( filename[2], READ_ONLY, TEXT_FILE)

	# Save item number if one was specified
	if ( nlist < 0 )
	   item = - nlist
	else
	   item = 0

	# Loop through entries
	nlist = 1
	jc = 1
	if (fd != 0) {
	   while ( getline(fd, linebuff) != EOF ) {
	      ic = 1
	      jc = 1
	      # Strip off leading white space
	      while (linebuff[jc] != EOS && jc <= SZ_LINE ) {
	         if (linebuff[jc] == ' ' || linebuff[jc] == '\t' ) 
	            jc = jc + 1
	         else
	            break
	      }

	      # Skip comments
	      if (linebuff[jc] == '#')
		 next

	      # No specific entry is desired so save each line
	      if ( item == 0 ) {

		 # Check to make sure list is not too long
		 if (nlist > maxlist) {
		    call sprintf(Memc[errmsg], SZ_LINE, toolong)
		    call pargstr( filename )
		    call pargi( maxlist )
		    call error( 1, Memc[errmsg] )
		 }

	         # Loop to EOS or end of line
	         while ( linebuff[jc] != EOS && jc <= SZ_LINE ) {

	            # If we hit a newline break and put EOS
	            if (linebuff[jc] == '\n' )
	               break

	            list[ic,nlist] = linebuff[jc]
	            ic = ic + 1
	            jc = jc + 1

	         }

	         # Put EOS and increment nlist if line is not blank
	         list[ic,nlist] = EOS
	         if ( strlen( list[1,nlist] ) > 0 )
	            nlist = nlist + 1

	      # Just want a single entry (the item-th)
	      } else if ( item > 0 ) {
	         if ( nlist < item ) {
	            nlist = nlist + 1
	            next
	         } else if ( nlist == item ) {
	            call strcpy( linebuff, list[1,1], SZ_LINE)
	            break
	         } else
	            break
	      }  # endif (item > 0)
	   }     # endwhile

	# Just a name not an @list
	} else {
	   call strcpy( filename, list[1,1], SZ_LINE )
	   nlist = 2	# ( we subtract 1 below )
	}

	# Decrement nlist by one
	nlist = nlist - 1
	
	# If we picked out a certain item, send back specification
	if (item != 0 )
	   nlist = item

	if (nlist == 0) {
	   call sprintf(Memc[errmsg], SZ_LINE, noentries)
	      call pargstr( filename )
	      call error( 1, Memc[errmsg] )
	}

	if (fd != 0)
	   call close(fd)
	call mfree(errmsg, TY_CHAR)
end
