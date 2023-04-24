# ft_fname -- append letter to file name
# The string RI is appended to the root portion of a file name.
# The file (image) name may include a logical directory, root name,
# extension, "cluster" index (group number), and/or image section.
# The name must not include a VMS directory specification because
# that would be confused with an image section.
#
# Phil Hodge, 20-Jul-1988  Subroutine created
# Phil Hodge, 24-Apr-1989  Change calling sequence of imparse.

procedure ft_fname (inname, ri, outname, maxch)

char	inname[ARB]		# i: input image name
char	ri[ARB]			# i: "r" or "i" to be appended to name
char	outname[ARB]		# o: name with "r" or "i" appended
int	maxch
#--
pointer sp
pointer cluster			# scratch for directory + root + extension
pointer ksection		# scratch (ignored)
pointer section			# scratch for image section (ignored)
pointer ldir, root		# scratch for directory, root name
int	cl_index, cl_size
int	dirlen, rootlen, flen
int	fnldir(), fnroot()

begin
	call smark (sp)
	call salloc (cluster, maxch, TY_CHAR)
	call salloc (ksection, maxch, TY_CHAR)
	call salloc (section, maxch, TY_CHAR)
	call salloc (ldir, maxch, TY_CHAR)
	call salloc (root, maxch, TY_CHAR)

	# Extract the "cluster" name, which is the logical directory,
	# file name root, and extension.  Ignore image section, etc.
	call imparse (inname,
		Memc[cluster], maxch, Memc[ksection], maxch,
		Memc[section], maxch, cl_index, cl_size)
	dirlen  = fnldir (Memc[cluster], Memc[ldir], maxch)
	rootlen = fnroot (Memc[cluster], Memc[root], maxch)

	# inname[1:flen] includes the logical directory and root but not
	# the (optional) dot and extension or cluster index or image section.
	flen = dirlen + rootlen
	call strcpy (inname, outname, flen)
	call strcat (ri, outname, maxch)
	call strcat (inname[flen+1], outname, maxch)

	call sfree (sp)
end
