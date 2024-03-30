include	"simtwo.h"

#* HISTORY *
#* B.Simon	14-Apr-95	original

# ZODIACAL -- Compute background counts from zodiacal light

procedure zodiacal (zodtab, ra, dec, jd, exptime, apscale, grftable, 
		    cmptable, nwave, wave, thruput, zspec, zcount)

char	zodtab[ARB]	# i: name of zodiacal light table
double	ra		# i: right ascension of target
double	dec		# i: declination of target
double	jd		# i: julian date of observation
real	exptime		# i: exposure time
double	apscale		# i: aperture scale
char	grftable[ARB]	# i: graph table name
char	cmptable[ARB]	# i: component lookup table name
int	nwave		# i: length of wavelength and thruput tables
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
real	zspec[ARB]	# o: zodiacal light spectrum
real	zcount		# o: counts / pix from zodiacal light
#--
int	degree
pointer	sp, expr, code, flux
real	pixlen, area, zflux

string	zband     ZOD_BAND
string	solar     "simulators$data/solar.dat"
string	expfmt    "%g * rn (%s, %s, 10.0, stmag)"

real	effstim()

begin
	# Check for null file

	if (zodtab[1] == EOS) {
	    zcount = 0.0
	    return
	}

	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (expr, SZ_LINE, TY_CHAR)
	call salloc (code, SZ_LINE, TY_INT)
	call salloc (flux, nwave, TY_REAL)

	# Compute pixel area in square degrees

	pixlen = 1.0
	call angtodegr (ZOD_UNITS, pixlen)
	pixlen = apscale / pixlen
	area = pixlen * pixlen

	# Compute flux in units of 10th magnitude suns

	call rdzodiac (zodtab, ra, dec, jd, zflux)
	zflux = zflux * exptime * area

	# Use expression evaluator to convert to more useful units

	call sprintf (Memc[expr], SZ_LINE, expfmt)
	call pargr (zflux)
	call pargstr (solar)
	call pargstr (zband)

	call syncompile (Memc[expr], Memi[code], SZ_LINE)
	call syncalc (Memi[code], SZ_LINE, NULL, nwave, wave, 
		      grftable, cmptable, Memr[flux], degree)

	# Compute counts from observation mode throughput and flux

	zcount = effstim (nwave, wave, thruput, Memr[flux], "counts")

	# Compute zodiacal light spectrum

	call amulr (thruput, Memr[flux], zspec, nwave)

	call sfree (sp)
end
