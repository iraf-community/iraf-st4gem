# All the routines in this file are obsolete.
# Last changed by Phil Hodge.

procedure load_wcstruct (im)

pointer	im	# i: image pointer
#--

begin
	# This subroutine presumably doesn't exist, so linking with
	# load_wcstruct should give an unresolved reference.
	call xxxxxx (im)
end

procedure get_wc_matrix (im, nax, wc_matrix)

pointer	im			# i: Image pointer
int	nax			# i: Number of axes in image
real	wc_matrix[nax,nax]	# i: Storage for the CD matrix
#--

begin
	call xxxxxx (im)
end

procedure reset_wcstruct (im)

pointer	im	# i: image pointer
#--

begin
	call xxxxxx (im)
end

procedure put_wcstruct( in, inaxis, out, outaxis)

pointer	in		# i: input  image pointer
int	inaxis		# i: input  image axis
pointer	out		# i: output image pointer
int	outaxis		# i: output image axis
#--

begin
	call xxxxxx (in)
end
