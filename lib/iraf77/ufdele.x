# ufdele -- delete a file
#
# Phil Hodge, 21-Feb-91  Subroutine created.

procedure ufdele (f77nam, istat)

				# i: file name
%      character*(*) f77nam
int	istat			# o: status return
#--
pointer sp
pointer fname			# scratch for file name
int	slen			# length of input string
int	errcode()

begin
	call smark (sp)
%      slen = len (f77nam)
	call salloc (fname, slen, TY_CHAR)

	call f77upk (f77nam, Memc[fname], slen)

	iferr {
	    call delete (Memc[fname])
	} then {
	    istat = errcode()
	} else {
	    istat = OK
	}

	call sfree (sp)
end
