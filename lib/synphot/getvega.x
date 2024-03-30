include	"libsynphot.h"

# GETVEGA -- Read a spectrum of Vega, sampled on the wavelength set

procedure getvega (nwave, wave, spec)

int     nwave           # i: length of wavelength and spectrum arrays
real    wave[ARB]       # i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectrum flux
#--

begin
	call rdstospec (VEGA, nwave, wave, spec)
end
