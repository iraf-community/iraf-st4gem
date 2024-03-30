#* HISTORY *
#* B.Simon	14-Apr-95	original
#* B.Simon	19-Jun-95	added earthlight and thermal contributions
#* B.Simon	30-Jun-95	spit into outbackgd and inbackgd

# INBACKGD -- Add background light from inside telescope to output image


procedure inbackgd (thermtab, exptime, apscale, nwave, wave, thruput, counts)
		    

char	thermtab[ARB]	# i: thermal background table
real	exptime		# i: exposure time
double	apscale		# i: aperture scale
int	nwave		# i: length of wavelength and thruput tables
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
real	counts		# i: total background counts
#--
real	tcount

begin
	# Thermal contribution to background

	call thermal (thermtab, exptime, apscale, 
		      nwave, wave, thruput, tcount)

	counts = tcount
end
