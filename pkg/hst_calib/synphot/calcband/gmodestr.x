# GMODESTR -- Get the rest of the mode string.  This routine parses the 
# command string until it hits a known compspec or compband command or
# an EOS.

procedure gmodestr( cmdstr, mode, iw, specdict, banddict )

char	cmdstr[ARB]	# i: String containg mode and other commands
char	mode[ARB]	# io: Parsed mode string.  First word is present
int	iw		# io: position in script
char	specdict[ARB]	# i: List of valid compspec commands
char	banddict[ARB]	# i: List of valid compband commands

#--

int	iw_old
int	ctowrd(), strdic()
char	word[SZ_LINE]

begin

	# Save position in script for pushback
	iw_old = iw

	# Loop while a non-zero length string is returned
	while( ctowrd( cmdstr, iw, word, SZ_LINE ) > 0 ) {

	   # If parsed word is a valid compspec or compband command then
	   # push it back and exit.
	   if ( strdic( word, word, SZ_LINE, specdict) > 0 ||
	        strdic( word, word, SZ_LINE, banddict) > 0 ) {
	      iw = iw_old
	      return

	   # If word is not a recognized command then concatinate it to the
	   # mode string, remember current position in cmdstr and continue.
	   } else {
	      call strcat( ",", mode, SZ_LINE )
	      call strcat( word, mode, SZ_LINE )
	      iw_old = iw
	   }
	}
end
