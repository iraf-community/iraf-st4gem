# BOXFUNC -- Box function

procedure boxfunc (center, fwhm, nwave, wave, band)

real	center		# i: center wavelength of box
real	fwhm		# i: box width
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength array
real	band[ARB]	# o: bandpass array
#--
int	iwave, iflag
real	shift[3]

string	badwidth  "Negative argument to function"
errchk	synphoterr

begin
	if (fwhm < 0.0)
	    call synphoterr (badwidth, "box")

	shift[1] = center + 0.5 * fwhm
	shift[2] = center - 0.5 * fwhm
	shift[3] = 0.0

	iflag = 0
	do iwave = nwave, 1, -1 {
	    if (wave[iwave] == shift[1])
		iflag = iflag + 1
	    if (wave[iwave] < shift[iflag+1])
		iflag = iflag + 1

	    if (iflag == 1) {
		band[iwave] = 1.0
	    } else {
		band[iwave] = 0.0
	    }
	}

end
