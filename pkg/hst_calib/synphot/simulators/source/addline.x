#* HISTORY *
#* B.Simon	18-Jul-95	original
#* B.Simon	15-Nov-95	better allocation of flux to pixels
#* B.Simon	12-Dec-95	modified for rectangular flux pixels

# ADDLINE -- Add a spectral line to the output image

procedure addline (xl, yl, grating, nwave, ypos, 
		   wave, line, nix, niy, out)

real	xl		# i: x center of line
real	yl		# i: y center of line
double	grating[ARB]	# i: grating dispersion coefficients
int	nwave		# i: length of spectral line
int	ypos[ARB]	# i: pixel location of spectrum
real	wave[ARB]	# i: wavelength set
real	line[ARB]	# i: spectral line
int	nix		# i: x dimension of output image
int	niy		# i: y dimension of output image
real	out[nix,niy]	# u: output image
#--
int	iwave, ix, jx, jy, kx
pointer	sp, lpos, rpos
real	xside, xpos, ldist, rdist, total, matrix[3]

string	badmatrix  "addline: inconsistent matrix calculation"
string	debugfmt   "jx = %d, jy = %d, line = %g, l = %f, c = %f, r = %f\n"

real	asumr()

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (lpos, nwave, TY_REAL)
	call salloc (rpos, nwave, TY_REAL)

	# Calculate cross dispersion of spectrum

	xside = xl - 0.5
	call crossdisp (xside, yl, grating, nwave, wave, Memr[lpos])

	xside = xl + 0.5
	call crossdisp (xside, yl, grating, nwave, wave, Memr[rpos])

	# Add spectrum to output image

	do iwave = 1, nwave {
	    xpos = 0.5 * (Memr[lpos+iwave-1] + Memr[rpos+iwave-1])

	    jx = xpos + 0.5
	    jy = ypos[iwave]

	    # Calculate what fraction of flux from pixel falls on the
	    # adjacent pixels

	    # Approximate formula assumes square flux area sides
	    # are aligned with pixel sides

	    ldist = (real (jx) - 0.5) - Memr[lpos+iwave-1]
	    rdist = Memr[rpos+iwave-1] - (real (jx) + 0.5)
	    
	    matrix[1] = max (ldist, 0.0)
	    matrix[3] = max (rdist, 0.0)
	    matrix[2] = 1.0 + min (ldist, 0.0) + min (rdist, 0.0)

	    total = asumr (matrix, 3)
	    call adivkr (matrix, total, matrix, 3)

	    do ix = 1, 3 {
		kx = jx + ix - 2
		if (kx < 1 || kx > nix)
		    next

		out[kx,jy] = out[kx,jy] + line[iwave] * matrix[ix]
	    }
	}

	call sfree (sp)
end

