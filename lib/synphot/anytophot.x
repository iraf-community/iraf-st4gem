include "libsynphot.h"

# ANYTOPHOT -- Convert flux units to photlam (photons cm-2 s-1 A-1)

int procedure anytophot (units, nwave, wave, flux)

char	units[ARB]	# i: input flux units
int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: wavelengths, in angstroms
real	flux[ARB]	# u: flux
#--
int	done, utype, iwave
pointer	sp, form, vflux
real	dwave, factor, area

string	badunits "Unknown flux units"
string	formlist FORMSTR

int	word_match()

begin
	call smark (sp)
	call salloc (form, SZ_FNAME, TY_CHAR)

	call strcpy (units, Memc[form], SZ_FNAME)
	call strfix (Memc[form])

	if (Memc[form] == EOS) {
	    done = NO

	} else {
	    done = YES
	    utype = word_match (Memc[form], formlist)

	    switch (utype) {
	    case 1:	# photlam
		;

	    case 2:	# counts
		call get_hstarea (area)

		dwave = abs (wave[2] - wave[1])
		flux[1] = flux[1] / (dwave * area)
		
		do iwave = 2, nwave-1 {
		    dwave = abs (wave[iwave+1] - wave[iwave-1]) / 2.0
		    flux[iwave] = flux[iwave] / (dwave * area)
		}

		dwave = abs (wave[nwave] - wave[nwave-1])
		flux[nwave] = flux[nwave] / (dwave * area)

	    case 3:	# flam
		factor = 1.0 / (H * C)
		do iwave = 1, nwave
		    flux[iwave] = factor * wave[iwave] * flux[iwave]

	    case 4:	# fnu
		factor = 1.0 / H
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] / wave[iwave]
	    case 5:	# photnu
		factor = C
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] / (wave[iwave] * 
							  wave[iwave])
	    case 6:	# jy
		factor = 1.0e-23 / H
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] / wave[iwave]

	    case 7:	# mjy
		factor = 1.0e-26 / H
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] / wave[iwave]

	    case 8:	# abmag
		factor = 1.0 / H
		do iwave = 1, nwave 
		    flux[iwave] = factor / wave[iwave] *
				  10.**( -0.4 * ( flux[iwave] - ABZERO ))

	    case 9:	# stmag
		factor = 1.0 / (H * C)
		do iwave = 1, nwave
		    flux[iwave] = factor * wave[iwave] *
				  10.**( -0.4 * (flux[iwave] - STZERO ))

	    case 10:	# vegamag
		call salloc (vflux, nwave, TY_REAL)
		call rdstospec (VEGA, nwave, wave, Memr[vflux])

		do iwave = 1, nwave
		    flux[iwave] = Memr[vflux+iwave-1] * 
				  10.0 ** (-0.4 * flux[iwave])

	    case 11:	# obmag
		call get_hstarea (area)

		dwave = abs (wave[2] - wave[1])
		flux[1] = 10.0 ** (-0.4 * flux[1]) / (dwave * area)
		
		do iwave = 2, nwave-1 {
		    dwave = abs (wave[iwave+1] - wave[iwave-1]) / 2.0
		    flux[iwave] = 10.0 ** (-0.4 * flux[iwave]) / (dwave * area)
		}

		dwave = abs (wave[nwave] - wave[nwave-1])
		flux[nwave] = 10.0 ** (-0.4 * flux[nwave]) / (dwave * area)

	    default:	# unknown units
		call synphoterr (badunits, units)

	    }
	}

	call sfree (sp)
	return (done)
end
