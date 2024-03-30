#* HISTORY *
#* B.Simon	03-Oct-94	original
#* B.Simon	01-Mar-00	added peaklam2

# PEAKLAM -- Find the peak thoughput of a passband

real procedure peaklam (nwave, filt)

int     nwave           # i: number of wavelengths
real    filt[ARB]       # i: throughput array
#--
int	iwave
real	peak

begin
	peak = filt[1]
	do iwave = 2, nwave {
	    if (filt[iwave] > peak)
		peak = filt[iwave]
	}

	return (peak)
end

procedure peaklam2 (nwave, wave, filt, peakwave, peakthru)

int     nwave           # i: number of wavelengths
real    wave[ARB]       # i: wavelength array
real    filt[ARB]       # i: throughput array
real	peakwave	# o: peak wavelength
real	peakthru	# o: peak throughput
#--
int	iwave, jwave

begin
	jwave = 1
	do iwave = 2, nwave {
	    if (filt[iwave] > filt[jwave])
		jwave = iwave
	}

	peakwave = wave[jwave]
	peakthru = filt[jwave]

end
