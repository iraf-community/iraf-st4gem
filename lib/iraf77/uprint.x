# uprint -- call printf
# This routine takes a Fortran character string, writes it to STDOUT,
# and flushes STDOUT.
#
# Phil Hodge,  6-Jun-1991  Subroutine created.

procedure uprint (f77str)

%      character*(*) f77str
#--
pointer sp
pointer buf		# scratch for SPP version of input string
int	length		# allocated length of Fortran string

begin
	# Add one for the \n that we will append.
%      length = len (f77str) + 1

	call smark (sp)
	call salloc (buf, length, TY_CHAR)
	call f77upk (f77str, Memc[buf], length)

	call strcat ("\n", Memc[buf], length)

	call printf (Memc[buf])
	call flush (STDOUT)

	call sfree (sp)
end
