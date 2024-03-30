include "source.h"
define	CHECKMASK		NO

#* HISTORY *
#* B.Simon	19-Jul-95	original

# SPECBACKGD -- Comput background for spectral model

procedure specbackgd (out, nix, niy, zodtab, ra, dec, jd, earthtab, eshine, 
		      exptime, apshape, apscale, grating, grftable, cmptable, 
		      nwave, wave, thruput)

real	out[nix,niy]	# u: output image buffer
int	nix		# i: first dimension of buffer
int	niy		# i: second dimension of buffer
char	zodtab[ARB]	# i: name of zodiacal light table
double	ra		# i: right ascension of target
double	dec		# i: declination of target
double	jd		# i: julian date of observation
char	earthtab[ARB]	# i: earthlight spectrum
real	eshine		# i: fraction of maximum earthlight
real	exptime		# i: exposure time
pointer	apshape		# i: aperture shape descriptor
double	apscale		# i: aperture scale
real	grating[ARB]	# i: grating dispersion parameters
char	grftable[ARB]	# i: graph table name
char	cmptable[ARB]	# i: component lookup table name
int	nwave		# i: length of wavelength and thruput tables
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
#--
int	ntx, nty
pointer	sp, zspec, espec, src, mask
real	zcount, ecount, counts, ox, oy

int	makemask()
real	asumr()

begin
	# Allocate memory for spectra

	call smark (sp)
	call salloc (zspec, nwave, TY_REAL)
	call salloc (espec, nwave, TY_REAL)

	# Allocate memory for source structure

	call salloc (src, SZ_SRCSTRUCT, TY_REAL)
	call aclri (Memi[src], SZ_SRCSTRUCT)

	call salloc (SRC_SCALARS(src), SZ_SRCSCALARS, TY_REAL)
	call salloc (SRC_WAVPTR(src), nwave, TY_REAL)
	call salloc (SRC_SPECPTR(src), nwave, TY_REAL)

	# Zodiacal light contribution to background

	call zodiacal (zodtab, ra, dec, jd, exptime, apscale, grftable, 
		       cmptable, nwave, wave, thruput, Memr[zspec], zcount)

	# Earth light contribution to background

	call earthshine (earthtab, eshine, exptime, apscale, 
			 nwave, wave, thruput, Memr[espec], ecount)

	counts = zcount + ecount

	# Fill in source structure

	ox = 0.5 * (nix + 1)
	oy = 0.5 * (niy + 1)
	SRC_XPOS(src) = ox
	SRC_YPOS(src) = oy

	SRC_NWAVE(src) = nwave
	call amovr (wave, SRC_WAVE(src), nwave)
	call aaddr (Memr[zspec], Memr[espec], SRC_SPEC(src), nwave)

	# Create aperture mask

	call malloc (mask, nix*niy, TY_REAL)

	if (makemask (apshape, apscale, ox, oy, nix, niy, 
		      nix, niy, Memr[mask]) == YES) {

	    # Diagnostic test of makemask

	    if (CHECKMASK == YES) {
		call amovr (Memr[mask], out, nix*niy)
		call sfree (sp)
		return
	    }

	    ntx = nix
	    nty = niy
	    call chopmask (ntx, nty, mask)

	    # Scale the mask by the background flux

	    call amulkr (Memr[mask], counts, Memr[mask], ntx*nty)
	    SRC_FLUX(src) = asumr (Memr[mask], ntx*nty)

	    # Add the background to the output

	    call addsource (src, grating, ntx, nty, nix, niy, Memr[mask], out)
	}

	call mfree (mask, TY_REAL)
	call sfree (sp)

end
