include <imhdr.h>
include	<imset.h>

# ft_map_s -- create a scratch image
# This routine creates an image as NEW_IMAGE with specified axis lengths
# up to 2-D, default (real) datatype.  A pseudo-random string will be
# appended to the prefix (SEED) to form the name.
#
# Phil Hodge, 12-Aug-1988  Subroutine created.

procedure ft_map_s (seed, nax, len1, len2, im)

char	seed[ARB]	# i: prefix for the image name
int	nax		# i: dimension of image; 1 or 2
int	len1		# i: length of first axis
int	len2		# i: length of second axis; ignored if nax=1
pointer im		# o: pointer to image header structure
#--
pointer sp
pointer imname		# scratch for image name
pointer immap()

begin
	call smark (sp)
	call salloc (imname, SZ_FNAME, TY_CHAR)

	call mktemp (seed, Memc[imname], SZ_FNAME)
	im = immap (Memc[imname], NEW_IMAGE, 0)
	IM_LEN(im,1) = len1
	if (nax > 1)
	    IM_LEN(im,2) = len2
	IM_PIXTYPE(im) = TY_REAL

	call sfree (sp)
end

# ft_reopen_s -- reopen a scratch image read/write
# This routine closes an image and reopens it read/write.  This allows
# the image to be used as scratch space.
#
# Phil Hodge, 12-Aug-1988  Subroutine created.

procedure ft_reopen_s (im)

pointer im		# io: pointer to image header structure
#--
pointer sp
pointer imname		# scratch for image name
pointer immap()

begin
	call smark (sp)
	call salloc (imname, SZ_FNAME, TY_CHAR)

	call imstats (im, IM_IMAGENAME, Memc[imname], SZ_FNAME)
	call imunmap (im)
	im = immap (Memc[imname], READ_WRITE, 0)

	call sfree (sp)
end

# ft_del_s -- delete a scratch image
# This routine deletes an image that was created for scratch space.
#
# Phil Hodge, 12-Aug-1988  Subroutine created.

procedure ft_del_s (im)

pointer im		# io: pointer to image header structure
#--
pointer sp
pointer imname		# scratch for image name

begin
	call smark (sp)
	call salloc (imname, SZ_FNAME, TY_CHAR)

	call imstats (im, IM_IMAGENAME, Memc[imname], SZ_FNAME)
	call imunmap (im)
	call imdelete (Memc[imname])

	call sfree (sp)
end
