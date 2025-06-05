# ufplin -- write a line to a text file
# This routine writes one line of text to a file.
# Note that a newline character is appended to the string before writing it.
#
# Phil Hodge, 21-Feb-91  Subroutine created.

procedure ufplin (fd, f77txt, istat)

int	fd			# i: file descriptor
				# i: the line of text
%      character*(*) f77txt
int	istat			# o: status return
#--
pointer sp
pointer cbuf			# scratch for SPP version of output string
int	slen			# length of output string
int	errcode()

begin
	call smark (sp)

	# One extra character for the '\n' which we will append to the string.
%      slen = len (f77txt) + 1
	call salloc (cbuf, slen, TY_CHAR)

	call f77upk (f77txt, Memc[cbuf], slen)
	call strcat ("\n", Memc[cbuf], slen)

	iferr {
	    call putline (fd, Memc[cbuf])
	} then {
	    istat = errcode()
	} else {
	    istat = OK
	}

	call sfree (sp)
end
