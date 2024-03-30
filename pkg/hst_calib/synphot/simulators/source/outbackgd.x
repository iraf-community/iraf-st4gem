#* HISTORY *
#* B.Simon	14-Apr-95	original
#* B.Simon	19-Jun-95	added earthlight contribution
#* B.Simon	30-Jun-95	spit into outbackgd and inbackgd

# OUTBACKGD -- Add background light from outside telescope to output image

procedure outbackgd (zodtab, ra, dec, jd, earthtab, eshine, exptime, apscale, 
		     grftable, cmptable, nwave, wave, thruput, counts)

char	zodtab[ARB]	# i: name of zodiacal light table
double	ra		# i: right ascension of target
double	dec		# i: declination of target
double	jd		# i: julian date of observation
char	earthtab[ARB]	# i: earthlight spectrum
real	eshine		# i: fraction of maximum earthlight
real	exptime		# i: exposure time
double	apscale		# i: aperture scale
char	grftable[ARB]	# i: graph table name
char	cmptable[ARB]	# i: component lookup table name
int	nwave		# i: length of wavelength and thruput tables
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
real	counts		# o: background counts
#--
pointer	sp, zspec, espec
real	zcount, ecount

begin
	# Allocate memory for spectra

	call smark (sp)
	call salloc (zspec, nwave, TY_REAL)
	call salloc (espec, nwave, TY_REAL)

	# Zodiacal light contribution to background

	call zodiacal (zodtab, ra, dec, jd, exptime, apscale, grftable, 
		       cmptable, nwave, wave, thruput, Memr[zspec], zcount)

	# Earth light contribution to background

	call earthshine (earthtab, eshine, exptime, apscale, 
			 nwave, wave, thruput, Memr[espec], ecount)

	counts = zcount + ecount
	call sfree (sp)
end
