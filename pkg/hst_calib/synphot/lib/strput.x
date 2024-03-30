# STRPUT -- Insert one string into another at a given index.  Truncate
# if to long

procedure strput( instr, mainstr, ic, maxch )

char	instr[ARB]	# i: inserted string
char	mainstr[ARB]	# io: string into which instr is inserted
int	ic		# i: position at which to insert
int	maxch		# i: length of mainstr	

int	jc
char	buff[SZ_LINE]

# Dave Bazell, Jan 1990

begin

	# Buffer mainstr from ic to end
	jc = 1	
	while( mainstr[ic+jc-1] != EOS && jc < SZ_LINE ) {
	   if ( ic+jc-1 < maxch ) {
	      buff[jc] = mainstr[ic+jc-1]
	      jc = jc + 1
	   } else
	      break
	}
	buff[jc] = EOS

	# Terminate mainstr at ic
	mainstr[ic] = EOS

	call strcat( instr, mainstr, maxch)
	call strcat( buff, mainstr, maxch)

end
