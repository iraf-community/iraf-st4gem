include <imio.h>
include <imhdr.h>
include "simtwo.h"

#* HISTORY *
#* B.Simon	21-Jun-95	original
#* B.Simon	30-Jun-95	derived from addbackfile

# FILEOPS -- Perform mathematical operation between input and output file

procedure fileops (out, nix, niy, infile, op)

real	out[nix,niy]	# u: output image buffer
int	nix		# i: first dimension of buffer
int	niy		# i: second dimension of buffer
char	infile[ARB]	# i: input file name
int	op		# i: mathematical operation
#--
int	iy
pointer	im, line

string	baddimen   "Input image is not two dimensional"
string	badsize    "Size of input image does not match detector"

int	checkdim()
pointer	immap(), imgl2r()

real    filezero() # Set result of illegal operation to zero
extern	filezero

begin
	# Check for input file

	if (infile[1] == EOS)
	    return

	# Open input image

	im = immap (infile, READ_ONLY, 0)

	# Check dimensions and size against output

	if (checkdim (im) != 2)
	    call printerr_str (baddimen, IM_NAME(im))

	if (IM_LEN(im,1) != nix || IM_LEN(im,2) != niy)
	    call printerr_str (badsize, IM_NAME(im))

	# Perform math operation on input image and output

	do iy = 1, niy {
	    line = imgl2r (im, iy)
	    switch (op) {
	    case ADDOP:
		call aaddr (out[1,iy], Memr[line], out[1,iy], nix)
	    case SUBOP:
		call asubr (out[1,iy], Memr[line], out[1,iy], nix)
	    case MULOP:
		call amulr (out[1,iy], Memr[line], out[1,iy], nix)
	    case DIVOP:
		call advzr (out[1,iy], Memr[line], out[1,iy], nix, filezero)
	    }
	}

	# Close input image

	call imunmap (im)
end

# FILEZERO -- Set result of illegal operation to zero

real procedure filezero (val)

real	val		# i: value operated on
#--

begin
	return (0.0)
end
