include	<iraf77.h>

# UPTLOG -- Put a message onto the log file

procedure uptlog (text, istat)

char 	text[ARB]
int	istat

char	logcmd[SZ_LINE]
char	loctext[SZ_LINE]
int	len, ip, op
int	strlen()

begin

    call strcpy ("putlog (\"", logcmd, SZ_LINE)

    # Escape any " in the message text
    len = strlen (text)
    op = 1
    for (ip=1; ip<=len; ip=ip+1) {
	if (text[ip] == '"') {
	    loctext[op] = '\\'
	    op = op + 1
	    loctext[op] = text[ip]
	    op = op + 1
	} else {
	    loctext[op] = text[ip]
	    op = op + 1
	}
    }
    loctext[op] = EOS

    # Put the message text into the logfile message
    call strcat (loctext, logcmd, SZ_LINE)
    call strcat ("\")", logcmd, SZ_LINE)

    # Send the command to the cl
    call clcmd (logcmd)
    istat = ER_OK
    return

end
