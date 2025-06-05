# ufclos -- close a file
# If the file is successfully closed, the value of fd will be set to zero.
#
# Phil Hodge, 21-Feb-91  Subroutine created.

procedure ufclos (fd, istat)

int	fd		# i: file descriptor
int	istat		# o: status return
#--
int	errcode()

begin
	iferr {
	    call close (fd)
	} then {
	    istat = errcode()
	} else {
	    fd = 0
	    istat = OK
	}
end
