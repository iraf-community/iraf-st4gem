include	<imhdr.h>

# ft_transpose -- transpose an image
# Divide the image into square blocks of size len_blk by len_blk.
# Transpose each block with an array transpose operator for single
# precision images.
#
# Phil Hodge, 30-Aug-1988  Subroutine copied from imtranspose.

procedure ft_transpose (im_in, im_out, len_blk)

pointer	im_in				# Input image descriptor
pointer	im_out				# Output image descriptor
int	len_blk				# 1D length of transpose block

int	x1, x2, nx
int	y1, y2, ny
pointer	buf_in, buf_out

pointer imgs2r(), imps2r()

begin
	# Break the input image into blocks of at most len_blk by len_blk.

	do x1 = 1, IM_LEN (im_in, 1), len_blk {
	    x2 = x1 + len_blk - 1
	    if (x2 > IM_LEN(im_in, 1))
	       x2 = IM_LEN(im_in, 1)
	    nx = x2 - x1 + 1

	    do y1 = 1, IM_LEN (im_in, 2), len_blk {
	        y2 = y1 + len_blk - 1
	        if (y2 > IM_LEN(im_in, 2))
		   y2 = IM_LEN(im_in, 2)
		ny = y2 - y1 + 1

		# Note that this is just for single-precision images.

		buf_in = imgs2r (im_in, x1, x2, y1, y2)
		buf_out = imps2r (im_out, y1, y2, x1, x2)
		call ft_tr2r (Memr[buf_in], Memr[buf_out], nx, ny)
	    }
	}
end


# ft_tr2r -- transpose a block of an image
# The arrays need not be identical.  (What does that mean?)

procedure ft_tr2r (a, b, nx, ny)

real	a[nx, ny], b[ny, nx]
int	nx, ny, x, y

begin
	do x = 1, nx
	   do y = 1, ny
	       b[y, x] = a[x, y]
end
