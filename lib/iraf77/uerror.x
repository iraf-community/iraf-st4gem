# UERROR -- call error
# This routine calls the iraf error subroutine.  The error number and
# error message are arbitrary except that the message will be truncated
# if it is longer than SZ_LINE.
#
# Phil Hodge, 23-May-1989  Subroutine created.

procedure uerror (errmes)

#				# i: text of error message
%      character*(*) errmes
#--
pointer sp
pointer message		# spp equivalent of error message
pointer omess		# buffer for original error message
int	iraf_code	# original error number
int	errget()

begin
	call smark (sp)
	call salloc (message, SZ_LINE, TY_CHAR)
	call salloc (omess, SZ_LINE, TY_CHAR)

	call f77upk (errmes, Memc[message], SZ_LINE)

	# Get the original error message and error code.  If there was no
	# error (i.e. if error was not called), iraf_code should be zero.
	iraf_code = errget (Memc[omess], SZ_LINE)
	if (iraf_code > 0) {
	    call strcat ("  ", Memc[message], SZ_LINE)
	    call strcat (Memc[omess], Memc[message], SZ_LINE)
	}
	call error (iraf_code, Memc[message])

	# There's no point in calling sfree.
end
