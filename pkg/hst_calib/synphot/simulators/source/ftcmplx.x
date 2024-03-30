# ft_cmplx -- Fourier transform of 2-D complex array
#
# Phil Hodge, 18-Jun-1992  Subroutine created.

procedure ft_cmplx (ximg, npix1, npix2, fwd)

complex ximg[npix1,npix2]	# io: array to be transformed
int	npix1			# i: length of first axis
int	npix2			# i: length of second axis
bool	fwd			# i: forward transform?
#--
pointer xwork			# scratch for complex array for a column
pointer trig			# scratch for array of cosines & sines
int	i, j			# loop indexes

begin
	# Allocate scratch space, and initialize the Fourier transform.
	call calloc (trig, 4*npix1 + 15, TY_REAL)
	call cffti (npix1, Memr[trig])

	# Do complex transform of each line in-place.
	if (npix1 > 1) {
	    if (fwd)				# forward transform
		do j = 1, npix2
		    call cfftf (npix1, ximg[1,j], Memr[trig])
	    else				# backward transform
		do j = 1, npix2
		    call cfftb (npix1, ximg[1,j], Memr[trig])
	}

	if (npix2 > 1) {

	    # If the axes are not the same length, we must reinitialize.
	    if (npix1 != npix2) {
		call mfree (trig, TY_REAL)
		call calloc (trig, 4*npix2 + 15, TY_REAL)
		call cffti (npix2, Memr[trig])
	    }
	    call malloc (xwork, npix2, TY_COMPLEX)

	    # Do complex transform of each column.
	    do i = 1, npix1 {			# do for each column

		do j = 1, npix2			# extract column into scratch
		    Memx[xwork+j-1] = ximg[i,j]

		if (fwd)			# forward transform
		    call cfftf (npix2, Memx[xwork], Memr[trig])
		else				# backward transform
		    call cfftb (npix2, Memx[xwork], Memr[trig])

		do j = 1, npix2			# copy scratch to column
		    ximg[i,j] = Memx[xwork+j-1]
	    }
	    call mfree (xwork, TY_COMPLEX)
	}
	call mfree (trig, TY_REAL)
end
