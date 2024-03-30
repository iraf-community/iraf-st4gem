#* HISTORY *
#* B.Simon	30-Mar-95	original
#* B.Simon	16-Jun-95	subtract mean from noise

# ADDNOISE -- Add poisson noise to the output image

procedure addnoise (out, nix, niy, quant, seed, noise)

real	out[nix,niy]	# u: output image buffer
int	nix		# i: first dimension of buffer
int	niy		# i: second dimension of buffer
bool	quant		# i: quantization error 
long	seed		# i: random number seed
int	noise[ARB]	# i: compiled noise expression
#--
int	iy
pointer	sp, mean, dev

string	noisemsg  "Adding noise to image"

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (mean, nix, TY_REAL)
	call salloc (dev, nix, TY_REAL)

	# Do computation one line at a time

	do iy = 1, niy {
	    # Print diagnostic message
	    call done_message (noisemsg, iy-1, niy)

	    # Evaluate noise expression to compute mean
	    call calcnoise (noise, nix, out[1,iy], Memr[mean])

	    # Compute deviant with this mean
	    call mknoise (seed, nix, Memr[mean], Memr[dev])

	    # Add deviant to output image
	    call aaddr (out[1,iy], Memr[dev], out[1,iy], nix)

	    # Subtract mean from noise
	    call asubr (out[1,iy], Memr[mean], out[1,iy], nix)

	    # Add quantization error

	    if (quant)
		call arndr (out[1,iy], out[1,iy], nix)
	}

	# Finish diagnostic message
	call done_message (noisemsg, niy, niy)

	call sfree (sp)
end
