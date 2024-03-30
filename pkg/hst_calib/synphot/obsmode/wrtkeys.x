# WRTKEYS -- Write the list of keywords to STDOUT

procedure wrtkeys (keyidx, keybuf, nkey, minlen)

int	keyidx[ARB]	# i: Array of keyword indices
char	keybuf[ARB]	# i: String buffer containing keywords
int	nkey		# i: Number of keywords
int	minlen		# i: Minimum allowable screen length
#--
int	ikey, pagelen, linelen

string	bl   " "
string	nl   "\n"

int	envgeti(), strlen()

begin
	# Return if nothing to print

	if (nkey == 0)
	    return

	# Sort keywords in alphabetical order

	call strsrt (keyidx, keybuf, nkey)

	# Get screen width

	pagelen = envgeti ("ttyncols")
	if (pagelen <= minlen)
	    pagelen = SZ_LINE

	linelen = 0
	do ikey = 1, nkey {
	    # Wrap keywords at end of line

	    linelen = strlen (keybuf[keyidx[ikey]]) + linelen + 1
	    if (linelen > pagelen) {
		linelen = strlen (keybuf[keyidx[ikey]]) + 1
		call putline (STDOUT, nl)
	    }

	    # Write keyword

	    call putline (STDOUT, keybuf[keyidx[ikey]])
	    call putline (STDOUT, bl)
	}

	call putline (STDOUT, nl)
end
