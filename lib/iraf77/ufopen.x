# ufopen -- open a text file
# This routine opens an ascii file.  Reasonable values of iomode include
# read-only (1), read-write (2), append (4), and new file (5).
#
# Phil Hodge, 21-Feb-91  Subroutine created.

procedure ufopen (f77nam, iomode, fd, istat)

				# i: file name
%      character*(*) f77nam
int	iomode			# i: I/O mode for file
int	fd			# o: file descriptor
int	istat			# o: status return
#--
pointer sp
pointer fname			# scratch for file name
int	slen			# length of input string
int	filetype		# file type
int	temp			# for the value of fd returned by open
int	open(), errcode()

begin
	call smark (sp)
%      slen = len (f77nam)
	call salloc (fname, slen, TY_CHAR)

	call f77upk (f77nam, Memc[fname], slen)

	filetype = TEXT_FILE
	iferr {
	    temp = open (Memc[fname], iomode, filetype)
	} then {
	    fd = 0
	    istat = errcode()
	} else {
	    fd = temp
	    istat = OK
	}

	call sfree (sp)
end
