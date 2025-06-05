include <fset.h>		# to check if output is redirected
include <iraf77.h>

define	PUT_TO_STDOUT	1
define	PUT_TO_STDERR	2
define	PUT_TO_LOGFILE	4

# UMSPUT -- Put a message to STDOUT or STDERR or user log file or any
# combination of those.  The message is limited to SZ_LINE in length.
# The message destination is coded as follows:
#
#  dest = 1:  standard output
#  dest = 2:  standard error
#  dest = 4:  log file, if keeplog=yes
#
# To put the message to more than one destination, add the codes.
#
#	?	?	   Original.
# Phil Hodge,  4-May-1994  Simplified.

procedure umsput (f77txt, dest, priority, istat)

% 	character*(*) f77txt	
int	dest			# i: destination
int	priority		# i: priority of message (ignored)
int	istat			# o: error return code, or OK
#--
char 	text[SZ_LINE]
bool	stdout_redir, stderr_redir
int	andi()
int	fstati()

begin
	istat = ER_OK

	# 1, 2, 4 (or the sum of those) are valid values for dest.
	if (dest < 0  || dest > 7)  {
	   istat = ER_UPTBADDEST
	   return
	}
	if (dest == MSG_NOP)		# MSG_NOP = 0
	   return

	# Convert character string to SPP string.
	call f77upk (f77txt, text, SZ_LINE)

	# Check where output has to go: 1 -> STDOUT, 2 -> STDERR,
	# 4 -> user logfile (if keeplog = yes).

	# Write to log file?
	if (andi (dest, PUT_TO_LOGFILE) > 0)
	    call uptlog (text, istat)

	call strcat ("\n", text, SZ_LINE)

	# Check if output to STDOUT is required.
	if (andi (dest, PUT_TO_STDOUT) > 0) {
	    call putline (STDOUT, text)
	    call flush (STDOUT)
	}

	# Check if output to STDERR is required.
	if (andi (dest, PUT_TO_STDERR) > 0) {

	    # Also write to STDOUT?  If so, we have already written the
	    # message, so if both STDOUT and STDERR are redirected, or
	    # if neither is redirected, we don't need to do it again.
	    if (andi (dest, PUT_TO_STDOUT) > 0) {
		stdout_redir = (fstati (STDOUT, F_REDIR) == YES)
		stderr_redir = (fstati (STDERR, F_REDIR) == YES)
		if (( stdout_redir &&  stderr_redir) ||
		    (!stdout_redir && !stderr_redir))
		    return
	    }

	    call putline (STDERR, text)
	}
end
