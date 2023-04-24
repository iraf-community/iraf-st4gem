# PIXORDER -- Comparison routine for mask pixels
#
# This procedure returns an integer indicating the order of two pixels, 
# according to the values stored in the order array. The ordering between
# pixels is determined by the precedence values of the pixels. The
# precedence value is the value the pixel would have had if bit flags were 
# sorted in precedence order. Precedence order means that the first bit flag 
# has the lowest precedence and the last bit flag has the highest precedence.
# The value is from this function follows this scheme:
#
#		if pix1 <  pix2, order <  0
#		if pix1 == pix2, order == 0
# 		if pix1 >  pix2, order >  0
#
# B.Simon	16-Jul-93	Original

int procedure pixorder (order, nflag, pix1, pix2)

int	order[ARB]	# i: Precedence order of flags
int	nflag		# i: Number of flags
int	pix1		# i: First pixel
int	pix2		# i: Second pixel
#--
int	ord, val1, val2

int	pixvalue()

begin
	if (nflag == 0) {
	    ord = pix1 - pix2

	} else {
	    if (pix1 == 0 && pix2 == 0) {
		ord = 0
	    } else if (pix1 == 0) {
		ord = -1
	    } else if (pix2 == 0) {
		ord = 1
	    } else {
		val1 = pixvalue (order, nflag, pix1)
		val2 = pixvalue (order, nflag, pix2)
		
		ord = val1 - val2
	    }
	}

	return (ord)
end

# PIXVALUE -- Get the ordered value of the pixel 

int procedure pixvalue (order, nflag, pix)

int	order[ARB]	# i: Precedence order of flags
int	nflag		# i: Number of flags
int	pix		# i: Pixelvalue
#--
int	iflag, value, bitval

begin
	# The ordered value of a pixel is the value the pixel 
	# would have if flags were in precedence order

	value = 0
	bitval = 1

	do iflag = 1, nflag {
	    if (pix < bitval)
		break

	    if (and (pix, bitval) != 0) {
		value = value + 2 ** (order[iflag] - 1)
	    }

	    bitval = 2 * bitval
	}

	return (value)
end
