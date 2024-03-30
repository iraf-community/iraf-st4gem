#* HISTORY *
#* B.Simon	03-Oct-94	original

# WIDTHLAM -- Calculate equivalent and rectangular width of a passband

procedure widthlam (nwave, wave, filt, equivwidth, rectwidth)

int     nwave           # i: number of wavelengths
real    wave[ARB]       # i: wavelength array
real    filt[ARB]       # i: throughput array
real	equivwidth	# o: equivalent width of passband
real	rectwidth	# o: width of rect with same peak & area as passband
#--
real	peak

real	peaklam(), syntegral()

begin
	# The peak throughput is the maximum value of the passband
	# The equivalent width of the passband is the integrated area 
	# under the passband divided bythe peak throughput

	peak = peaklam (nwave, filt)
	equivwidth = syntegral (nwave, wave, filt)

	if (peak > 0.0) {
	    rectwidth = equivwidth / peak
	} else {
	    rectwidth = 0.0
	}

end
