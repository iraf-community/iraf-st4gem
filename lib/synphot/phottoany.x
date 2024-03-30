include "libsynphot.h"
define	SZ_VAL		31

#* HISTORY *
#* B.Simon	26-Jan-96	Better error message for vegamag
#* V.Laidler    14-Jan-05	Warn user when truncating negative fluxes

# PHOTTOANY -- Convert flux units from photlam (photons cm-2 s-1 A-1)

int procedure phottoany (units, nwave, wave, flux)

char	units[ARB]	# i: input flux units
int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: wavelengths, in angstroms
real	flux[ARB]	# u: flux
#--
int	done, utype, iwave, negcount
pointer	sp, form, value, vflux
real	dwave, factor, area

string	badunits  "Unknown flux units"
string	vegarange "Wavelength outside defined range of Vega's spectrum"
string	formlist  FORMSTR
string  negflux1 "Input spectrum contains flux values <= 0"
char	negmsg[SZ_LINE]


int	word_match()

begin
	call smark (sp)
	call salloc (form, SZ_FNAME, TY_CHAR)
	call salloc (value, SZ_VAL, TY_CHAR)

	call strcpy (units, Memc[form], SZ_FNAME)
	call strfix (Memc[form])

	negcount = 0

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
		flux[1] = flux[1] * dwave * area
		
		do iwave = 2, nwave-1 {
		    dwave = abs (wave[iwave+1] - wave[iwave-1]) / 2.0
		    flux[iwave] = flux[iwave] * dwave * area
		}

		dwave = abs (wave[nwave] - wave[nwave-1])
		flux[nwave] = flux[nwave] * dwave * area

	    case 3:	# flam
		factor = H * C
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] / wave[iwave]

	    case 4:	# fnu
		factor = H
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] * wave[iwave]

	    case 5:	# photnu
		factor = 1.0 / C
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] * 
				  wave[iwave] * wave[iwave]

	    case 6:	# jy
		factor = 1.0e23 * H
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] * wave[iwave]

	    case 7:	# mjy
		factor = 1.0e26 * H
		do iwave = 1, nwave
		    flux[iwave] = factor * flux[iwave] * wave[iwave]

	    case 8:	# abmag
		factor = H
		do iwave = 1, nwave {
		    flux[iwave] = factor * flux[iwave] * wave[iwave]

		    if (flux[iwave] <= 0.0) {
			flux[iwave] = 100.0
			negcount = negcount+1
		    } else {
			flux[iwave] = -2.5 * alog10 (flux[iwave]) + ABZERO
		    }
		}

	    case 9:	# stmag
		factor = H * C
		do iwave = 1, nwave {
		    flux[iwave] = factor * flux[iwave] / wave[iwave]

		    if (flux[iwave] <= 0.0) {
			flux[iwave] = 100.0
			negcount = negcount + 1
		    } else {
			flux[iwave] = -2.5 * alog10 (flux[iwave]) + STZERO
		    }
		}

	    case 10:	# vegamag

		call salloc (vflux, nwave, TY_REAL)
		call rdstospec (VEGA, nwave, wave, Memr[vflux])

		do iwave = 1, nwave {
		    if (flux[iwave] <= 0.0) {
			flux[iwave] = 100.0
			negcount = negcount + 1
		    } else if (Memr[vflux+iwave-1] <= 0.0){
			call sprintf (Memc[value], SZ_VAL, "%f")
			call pargr (wave[iwave])
			call synphoterr (vegarange, Memc[value])
		    } else {
			flux[iwave] = -2.5 * alog10 (flux[iwave] / 
						     Memr[vflux+iwave-1])
		    }
		}

	    case 11:	# obmag

		call get_hstarea (area)

		do iwave = 1, nwave {
		    if (iwave == 1) {
			dwave = abs (wave[2] - wave[1])
		    } else if (iwave == nwave) {
			dwave = abs (wave[nwave] - wave[nwave-1])
		    } else {
			dwave = abs (wave[iwave+1] - wave[iwave-1]) / 2.0
		    }

		    if (flux[iwave] <= 0.0) {
			flux[iwave] = 100.0
			negcount = negcount + 1
		    } else {
			flux[iwave] = -2.5 * alog10(flux[iwave] * dwave * area)
		    }
		}

	    default:	# unknown units
		call synphoterr (badunits, units)

	    }
	}

	if (negcount != 0) {

	  call sprintf(negmsg, SZ_LINE, " mag value set to 100 for %d of %d bins")
	  call pargi(negcount)
	  call pargi(nwave)

	  call synphotwarn(negflux1, negmsg) }

	call sfree (sp)
	return (done)
end
