# ufglin -- get a line from a text file
# This routine gets one line of text from a file.  Getlline (get long line)
# is used rather than getline so lines longer than SZ_LINE can be read.
# This routine (and getline) return a newline character at the end of the
# string, but we remove this newline before packing to a Fortran character
# string.  When the end of file is read, the output buffer is set to blank,
# and istat is set to EOF (-2).
#
# Phil Hodge, 21-Feb-91  Subroutine created.

procedure ufglin (fd, f77txt, istat)

int	fd			# i: file descriptor
				# o: the line of text
%      character*(*) f77txt
int	istat			# o: status return; zero is OK
#--
pointer sp
pointer cbuf			# scratch for SPP version of output string
int	slen			# allocated length of output string
int	extlen			# slen extended by SZ_LINE
int	status			# status return from getline or getlline
int	i			# loop index
int	errcode(), getlline()

begin
	call smark (sp)
%      slen = len (f77txt)

	# Add extra space to slen because the third argument of getlline
	# needs to be larger than the actual maxch by SZ_LINE.
	extlen = slen + SZ_LINE
	call salloc (cbuf, extlen, TY_CHAR)

	iferr {
	    status = getlline (fd, Memc[cbuf], extlen)

	} then {
	    istat = errcode()

	} else {
	    if (status == EOF) {
		istat = EOF
		Memc[cbuf] = EOS
	    } else {
		istat = OK
	    }

	    # Remove the '\n'.
	    do i = 0, extlen-1 {
		if (Memc[cbuf+i] == EOS)
		    break
		else if (Memc[cbuf+i] == '\n')
		    Memc[cbuf+i] = EOS
	    }
	}

	call f77pak (Memc[cbuf], f77txt, slen)

	call sfree (sp)
end
