include "../fourier.h"

# This file contains the following subroutines for allocating or
# releasing memory (in case the corresponding image does not exist)
# or getting/putting a line from/to an image (if the image does exist):
#
# ft_malloc
# ft_mfree
# ft_gl1r, ft_pl1r
# ft_gl2r, ft_pl2r
#
# Phil Hodge, 8-Aug-1988  Subroutine created.

# ft_malloc -- allocate memory
# If there is no input or output image for the real part or imaginary
# part, allocate memory to be used as if there were an image.  The
# memory should be set to zero if it's an input image.

procedure ft_malloc (ft, npts, clear_mem, Re_ptr, Im_ptr)

pointer ft		# i: FT structure
int	npts		# i: size of array(s) to allocate
int	clear_mem	# i: calloc rather than malloc if clear_mem=YES
pointer Re_ptr		# o: ptr to memory for real part, or unchanged
pointer Im_ptr		# o: ptr to memory for imaginary part, or unchanged
#--

begin
	if (FT_REPT(ft) == NULL) {	# is there no image for real part?
	    if (clear_mem == YES)
		call calloc (Re_ptr, npts, TY_REAL)
	    else
		call malloc (Re_ptr, npts, TY_REAL)
	}

	if (FT_IMPT(ft) == NULL) {	# is there no image for imaginary part?
	    if (clear_mem == YES)
		call calloc (Im_ptr, npts, TY_REAL)
	    else
		call malloc (Im_ptr, npts, TY_REAL)
	}
end

# ft_mfree -- free memory
# If the real (imaginary) part does not exist, free memory that was
# used in place of an input or output file.

procedure ft_mfree (ft, Re_ptr, Im_ptr)

pointer ft		# i: FT structure
pointer Re_ptr		# io: ptr to memory for real part
pointer Im_ptr		# io: ptr to memory for imaginary part
#--

begin
	if (FT_REPT(ft) == NULL)	# is there no image for real part?
	    call mfree (Re_ptr, TY_REAL)

	if (FT_IMPT(ft) == NULL)	# is there no image for imaginary part?
	    call mfree (Im_ptr, TY_REAL)
end

# ft_gl1r -- call imgl1r
# If the real (imaginary) 1-D input image does exist, get a pointer to it.

procedure ft_gl1r (realpt, imagpt, Re_ptr, Im_ptr)

pointer realpt		# i: ptr to imio structure for real part
pointer imagpt		# i: ptr to imio structure for imaginary part
pointer Re_ptr		# o: ptr to input image for real part
pointer Im_ptr		# o: ptr to input image for imaginary part
#--
pointer imgl1r()

begin
	if (realpt != NULL)
	    Re_ptr = imgl1r (realpt)

	if (imagpt != NULL)
	    Im_ptr = imgl1r (imagpt)
end

# ft_pl1r -- call impl1r
# If the real (imaginary) 1-D output image does exist, get a pointer to it.

procedure ft_pl1r (realpt, imagpt, Re_ptr, Im_ptr)

pointer realpt		# i: ptr to imio structure for real part
pointer imagpt		# i: ptr to imio structure for imaginary part
pointer Re_ptr		# o: ptr to output image for real part
pointer Im_ptr		# o: ptr to output image for imaginary part
#--
pointer impl1r()

begin
	if (realpt != NULL)
	    Re_ptr = impl1r (realpt)

	if (imagpt != NULL)
	    Im_ptr = impl1r (imagpt)
end

# ft_gl2r -- call imgl2r
# If the real (imaginary) 2-D input image does exist, get a pointer to
# one line of it.

procedure ft_gl2r (realpt, imagpt, line, Re_ptr, Im_ptr)

pointer realpt		# i: ptr to imio structure for real part
pointer imagpt		# i: ptr to imio structure for imaginary part
int	line		# i: which line of image to read
pointer Re_ptr		# o: ptr to input image for real part
pointer Im_ptr		# o: ptr to input image for imaginary part
#--
pointer imgl2r()

begin
	if (realpt != NULL)
	    Re_ptr = imgl2r (realpt, line)

	if (imagpt != NULL)
	    Im_ptr = imgl2r (imagpt, line)
end

# ft_pl2r -- call impl2r
# If the real (imaginary) 2-D output image does exist, get a pointer to
# one line of it.

procedure ft_pl2r (realpt, imagpt, line, Re_ptr, Im_ptr)

pointer realpt		# i: ptr to imio structure for real part
pointer imagpt		# i: ptr to imio structure for imaginary part
int	line		# i: which line of image to write
pointer Re_ptr		# o: ptr to output image for real part
pointer Im_ptr		# o: ptr to output image for imaginary part
#--
pointer impl2r()

begin
	if (realpt != NULL)
	    Re_ptr = impl2r (realpt, line)

	if (imagpt != NULL)
	    Im_ptr = impl2r (imagpt, line)
end
