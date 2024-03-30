# SELECTC -- select from a list of commands
#
# Interactive selection of a name from a fixed list
#
# Input:
#	PROMPT	= Character string to be used in the prompt
#	NLIST	= Number of names in the list
#	NAME	= List of names from which a selection is to be made
# Output:
#	NSELECT	= The sequence number of the selected name
#		= 0 if no name was selected.
#
# Oct 1985 by Keith Horne at STScI
# Mar 1988 KDH @ STScI - minor modification
# Jan 1989 E. Medeiros - SPP version
#
include <tbset.h>
procedure selectc ( prompt, nlist, name, nselect )

int	nmatch				# command group number matched
int	match[25]			# command group
int 	nlist				# number of commands in list
int	nselect				# sequence number of selected command
int	status				# STDIN read status
int	scan				# STDIN read function value
int	strmatch			# string match function value
int     stridx				# character match function value
int	i				# loop counter

char	name[SZ_LINE,ARB]		# command strings
char	prompt[SZ_LINE]			# command string list header
char	command[SZ_LINE]		# interactive response string

begin

	repeat { # until user request
	   # initialize the selection number
	   nselect = 0

	   # get the command token and convert command to upper case
	   call printf ( "%20s\n" )
	   call pargstr ( prompt )
	   call flush ( STDOUT )
	   status = scan()
	   call gargstr ( command, SZ_LINE )
	   call strupr ( command[1] )

	   # scan command name list for match
	   nmatch = 0
	   do i = 1, nlist { # sift throught command list for match
	      if ( stridx ( command[1], name[1,i] ) == 1 ) {
                 nmatch = nmatch + 1
                 match[nmatch] = i
              }
           }
	   if ( strmatch ( "QUIT", command ) > 1 ) { # accept state
	      nselect = 100
	      return
	   }
       	   if ( nmatch > 1 ) { #abiguous selections have been made
 	      call printf ( "Ambiguous selection\n" )

	   }else if ( nmatch <= 0 ) { # no match found interrogate user again
	      call printf ( "Valid commands\n" )
 	      do i = 1, nlist { 
	         call printf ( "%s\n" )
	         call pargstr ( name[1,i] )
	      }
	      call printf ( "OUIT ... Quit.\n" )

	   }else{ # valid command found
	      nselect = match[1]
	      return
	   }
	   call flush ( STDOUT )
        }
end
