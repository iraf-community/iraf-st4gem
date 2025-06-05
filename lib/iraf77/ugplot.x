include <iraf77.h>

# UGPLOT -- Plot a vector.  This routine is provided for the convenience of
# the user who does not need to exercise fine control over the details of
# how the plot is generated.

procedure ugplot (v, npts, x1, x2, f77ttl, istat)

real	v[ARB]			# data vector
int	npts			# number of data points
real	x1, x2			# range of X in data vector
% 	character*(*) f77ttl
int	istat

int	txtlen
pointer	sp, title

begin
	istat = ER_OK

	# Convert plot title to SPP name
	txtlen = len (f77ttl)
	call smark (sp)
	call salloc (title, txtlen, TY_CHAR)
	call f77upk (f77ttl, Memc[title], txtlen)

	iferr (call gplotv (v, npts, x1, x2, Memc[title]))
	   istat = ER_GRAPHPLOT

	call sfree (sp)
end
