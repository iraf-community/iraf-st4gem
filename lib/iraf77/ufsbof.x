# ufsbof -- seek the beginning of a file
# This routine rewinds a file opened read_only or read_write.  It will
# give unpredictable results if used on a newly created output file.
#
# Phil Hodge, 26-Feb-91  Subroutine created.

procedure ufsbof (fd, istat)

int	fd		# i: file descriptor
int	istat		# o: status return
#--
long	offset		# the location to seek
int	errcode()

begin
	offset = BOF

	iferr {
	    call seek (fd, offset)
	} then {
	    istat = errcode()
	} else {
	    istat = OK
	}
end
