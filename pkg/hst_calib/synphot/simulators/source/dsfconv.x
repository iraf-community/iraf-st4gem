include	"otf.h"

#* HISTORY *
#* B.Simon	11-Aug-95	based on fconvolve

# DSFCONV -- Convolve the output image with the detector spread function

procedure dsfconv (dsf, nix, niy, out)

pointer	dsf		# i: detector point spread function descriptor
int	nix		# i: first dimension of buffer
int	niy		# i: second dimension of buffer
real	out[nix,niy]	# u: output image buffer
#--
bool	fwd
complex	xnpix
int	npx, npy
pointer	buf1, buf2, adsf

begin
	# Check for null detector spread function

	if (dsf == NULL)
	    return

	# Allocate complex buffers for fourier transforms

	call malloc (buf1, nix*niy, TY_COMPLEX)
	call malloc (buf2, nix*niy, TY_COMPLEX)

	# Copy data into complex buffers

	call achtrx (out, Memx[buf1], nix*niy)

	npx = OTF_NXPIX(dsf)
	npy = OTF_NYPIX(dsf)
	adsf = OTF_BUFFER(dsf)
	call dsfshift (npx, npy, Memr[adsf], nix, niy, Memx[buf2])

	# Do the forward Fourier transform in-place in the complex arrays.

	fwd = true
	call ft_cmplx (Memx[buf1], nix, niy, fwd)
	call ft_cmplx (Memx[buf2], nix, niy, fwd)

	# Multiply the complex arrays, writing back to buf1.

	call amulx (Memx[buf1], Memx[buf2], Memx[buf1], nix*niy)

	# Normalize, working in-place in buf1.

	xnpix = complex (nix*niy, 0)	# both arguments are integer
	call adivkx (Memx[buf1], xnpix, Memx[buf1], nix*niy)

	# Take the inverse Fourier transform in-place in buf1.

	fwd = false
	call ft_cmplx (Memx[buf1], nix, niy, fwd)
		       
	# Copy complex buffer back into output array

	call achtxr (Memx[buf1], out, nix*niy)

	# Free buffers
       
	call mfree (buf1, TY_COMPLEX)
	call mfree (buf2, TY_COMPLEX)
end

