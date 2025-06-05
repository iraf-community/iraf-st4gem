include	<ctype.h>

# UERGET -- Get error message text associated with an error number.
#
# Nelson Zarate            subroutine created
# Phil Hodge, 1-June-1989  modified to also return the original error message
# Lucy Willard, 2-June-1989  comment out declaration of nchars and junk

procedure uerget (errcode, text)

int	errcode			# i: error message number
#				# o: error message
%	character*(*)	text

char	outstr[SZ_LINE]		# output string
char	errmsg[SZ_LINE]
char	buf[SZ_LINE]
#int	codelen, nchars, ip, junk
int	codelen, ip
int	k
int	iraf_code		# iraf error code
int	fd			# fd for iraf77errs file
int	strncmp(), itoc()
int	open()
int	getline()
int	errget()
define	nofile_ 91
define  noer_   92

begin
	# Get the original error message (buf) and error number (iraf_code).
	iraf_code = errget (buf, SZ_LINE)

	# Write the error number and message into the buffer outstr.
	# If the iraf error routine was not called, iraf_code should be
	# zero, and there will be no message.
	if (iraf_code == 0) {
	    outstr[1] = EOS
	} else {
	    call sprintf (outstr, SZ_LINE, "%d %s")
		call pargi (iraf_code)
		call pargstr (buf)
	}

	# Write error code to errmsg.
	codelen  = itoc (errcode, errmsg, SZ_LINE)

	# Open the error message file.
	iferr {
	    fd = open ("stsdas$lib/iraf77errs", READ_ONLY, TEXT_FILE)
	} then {
	    goto nofile_
	}

	# Read the error message file, and compare the number at the
	# beginning of the line with the input errcode.

	while (getline (fd, buf) != EOF) {
	    do k = 1, SZ_LINE {
		if (buf[k] == '\n')
		    buf[k] = EOS
	    }
	    if (strncmp (buf, errmsg, codelen) == 0) {
		# Skip the error code prefix and the blank which follows.
		for (ip=codelen+1;  IS_WHITE(buf[ip]);  ip=ip+1)
		    ;

		# Output stsdas system error message.
		if (outstr[1] != EOS)
		    call strcat ("; ", outstr, SZ_LINE)
		call strcat (buf[ip], outstr, SZ_LINE)
		break
	    }
	}
	call close (fd)

noer_
	call f77pak (outstr, text, SZ_LINE)
	return

nofile_
	call strcat (" (no stsdas error file)", outstr, SZ_LINE)
	call f77pak (outstr, text, SZ_LINE)
end
