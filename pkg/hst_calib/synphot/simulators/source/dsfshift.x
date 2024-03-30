#* HISTORY *
#* B.Simon	14-Aug-95	original

# DSFSHIFT -- Shift and zero pad a PSF for the fourier transform

procedure dsfshift (npx, npy, adsf, nix, niy, xdsf)

int	npx		# i: x dimension of detector PSF
int	npy		# i: y dimension of detector PSF
real	adsf[npx,npy]	# i: detector psf
int	nix		# i: x dimension of output PSF
int	niy		# i: y dimension of output PSF
complex	xdsf[nix,niy]	# o: output PSF, shifted and zero padded
#--
int	shiftx, shifty, iix, iiy, jix, jiy

begin
	# Compute shift needed to place center of psf in corner 

	shiftx = 1 - (npx + 1) / 2
	shifty = 1 - (npy + 1) / 2

	# Copy psf to output array

	do iiy = 1, niy {
	    jiy = iiy + shifty
	    if (jiy > niy)
		jiy = jiy - niy
	    if (jiy <= 0)
		jiy = jiy + niy

	    # Zero pad psf

	    call aclrx (xdsf[1,jiy], nix)

	    if (iiy <= npy) {
		do iix = 1, npx {
		    jix = iix + shiftx
		    if (jix > nix)
			jix = jix - nix
		    if (jix <= 0)
			jix = jix + nix

		    xdsf[jix,jiy] = adsf[iix,iiy]
		}
	    }
	}

end
