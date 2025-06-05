include	<ctype.h>
include <iraf77.h>

# SDASPARSE -- Parse an image specification into the cluster root, cluster 
# index, and cluster size,
#
#	Syntax:		cluster[cl_index/cl_size] or cluster[*]
#
# where all fields are optional except the cluster name.  In the limiting case
# (cl_size = 1) the cluster name and image name are the same.  CL_INDEX and
# CL_SIZE must be simple positive decimal integer constants, if given.  The
# [ character must be escaped to be included in the filename of the cluster.
# [cl_index/cl_size] may be replaced by [*] which means all images in a 
# multi-image file.
#
# NOTE: Throw away code -- Will not be needed when full blown IRAF template
#       expander will be available.
#       Based on imio$imparse.

procedure sdas_parse (imspec,
	cluster, cl_index, cl_size, istat)

char	imspec[SZ_PATHNAME]	# full image specification
char	cluster[SZ_PATHNAME]	# receives cluster name
int	cl_index		# receives cluster index (default 0)
int	cl_size			# receives cluster size (default 0)
int	istat			# return status
int	ch, n
int	ip, op, lbrack

begin
	istat = ER_OK
	ip = 1
	op = 1

	# Extract cluster name.  The first (unescaped) [ marks the start of
	# the cl_index subscript.

	for (ch=imspec[ip];  ch != EOS && ch != '[';  ch=imspec[ip]) {
	    if (ch == '\\' && imspec[ip+1] == '[') {
		cluster[op] = '\\'
		op = op + 1
		cluster[op] = '['
		ip = ip + 1
	    } else
		cluster[op] = ch

	    op = min (SZ_PATHNAME, op + 1)
	    ip = ip + 1
	}

	cluster[op] = EOS
	lbrack      = ip
	cl_index    = 0
	cl_size     = 0

	if (ch == EOS)
	    return

	# If we have a [...] field, find out if it is a wildcard specif.
	# A cl_index subscript is anything with the syntax
	# [ddd] or [ddd/ddd] or [*]

	ip = ip + 1
	n  = 0

	for (ch=imspec[ip];  ch != EOS && ch != ']';  ch=imspec[ip]) {
	    if (IS_DIGIT(ch)) {
		n = (n * 10) + TO_INTEG(ch)
	    } else if (ch == '/') {
		cl_index = max (n, 1)
		n = 0
	    } else if (ch == '*') {
		cl_size = -1
		return
	    } else {
		# bad image specification - incorrect syntax
		istat = ER_IMBADSEC
		return
	    }
	    ip = ip + 1
	}

	if (cl_index == 0)
	    cl_index = n
	else
	    cl_size = n

end
