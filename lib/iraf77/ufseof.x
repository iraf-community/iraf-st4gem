# ufseof -- seek the end of a file
# This routine goes to the end of a file opened read_only or read_write.
# There is no need to call it for a newly created file (new_file).
#
# Phil Hodge, 26-Feb-91  Subroutine created.

procedure ufseof (fd, istat)

int	fd		# i: file descriptor
int	istat		# o: status return
#--
long	offset		# the location to seek
int	errcode()

begin
	offset = EOF

	iferr {
	    call seek (fd, offset)
	} then {
	    istat = errcode()
	} else {
	    istat = OK
	}
end
