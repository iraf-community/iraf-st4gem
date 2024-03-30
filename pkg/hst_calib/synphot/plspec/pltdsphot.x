include	<gset.h>
include "dsphot.h"

define	IS_NONE		1
define	IS_POINT	2
define	IS_CONT		3
define	IS_HIST		4

define	HAS_HORIZ	1
define	HAS_VERT	2

define	ALL_BAR		3
define	NO_HORIZ	2
define	NO_VERT		1
define	NO_BAR		0

#* HISTORY *
#* B.Simon	13-Jun-94	original

# PLTDSPHOT -- Plot the spectrophotometry file

procedure pltdsphot (gp, dsphot, errtyp)

pointer	gp		# i: graphics descriptor
pointer	dsphot		# i: spectrophotometry descriptor
char	errtyp[ARB]	# i: type of error plot
#--
int	ic, nwave, line, bar, tbar, iflux
int	barscreen[4]
pointer	sp, wave, flux, err, fwhm, tflux, terror, tfwhm

string	linelist  "npcb"
string	barlist   "hv"

data	barscreen / NO_BAR, ALL_BAR, NO_VERT, NO_HORIZ /
string	badline   "Ambiguous choice of error type"

int	stridx()

begin
	if (dsphot == NULL)
	    return

	nwave = SPT_NWAVE(dsphot)
	wave = SPT_WAVE(dsphot)
	flux = SPT_FLUX(dsphot)
	err = SPT_ERROR(dsphot)
	fwhm = SPT_FWHM(dsphot)

	# Allocate memory for temporary arays

	call smark (sp)
	call salloc (tflux, nwave, TY_REAL)
	call salloc (terror, nwave, TY_REAL)
	call salloc (tfwhm, nwave, TY_REAL)

	# Get error line and bar type from error type string

	line = 0
	for (ic = 1; errtyp[ic] != EOS; ic = ic + 1)
	    line = line + stridx (errtyp[ic], linelist)

	if (line == 0)
	    line = IS_NONE

	if (line > IS_HIST)
	    call printerr_str (badline, errtyp)


	bar = 0
	for (ic = 1; errtyp[ic] != EOS; ic = ic + 1)
	    bar = bar + stridx (errtyp[ic], barlist)

	bar = and (bar, barscreen[line])

	do iflux = 1, SPT_NFLUX(dsphot) {
	    # Set temporary bar type according to existence of errors

	    tbar = bar
	    if (IS_INDEFR(Memr[err])) {
		tbar = and (tbar, NO_VERT)
	    } else {
		call amulkr (Memr[err], -1.0, Memr[terror], nwave)
	    }

	    if (IS_INDEFR(Memr[fwhm])) {
		tbar = and (tbar, NO_HORIZ)
	    } else {
		call amulkr (Memr[fwhm], -1.0, Memr[tfwhm], nwave)
	    }
		
	    # Plot flux and associated error flux, if any

	    switch (line) {
	    case IS_NONE:
		call gpline (gp, Memr[wave], Memr[flux], nwave)

	    case IS_POINT:
		if (and (tbar, HAS_HORIZ) > 0)
		    call pltbar (gp, Memr[wave], Memr[flux], nwave, GM_HLINE, 
				 Memr[tfwhm])

		if (and (tbar, HAS_VERT) > 0)
		    call pltbar (gp, Memr[wave], Memr[flux], nwave, GM_VLINE,
				  Memr[terror])

		if (tbar == NO_BAR)
		    call gpmark (gp, Memr[wave], Memr[flux], nwave, GM_POINT,
				 1.0, 1.0)

	    case IS_CONT:
		call gpline (gp, Memr[wave], Memr[flux], nwave)

		call aaddr (Memr[flux], Memr[err], Memr[tflux], nwave)
		call gpline (gp, Memr[wave], Memr[tflux], nwave)

		call asubr (Memr[flux], Memr[err], Memr[tflux], nwave)
		call gpline (gp, Memr[wave], Memr[tflux], nwave)

		if (and (tbar, HAS_HORIZ) > 0)
		    call pltbar (gp, Memr[wave], Memr[flux], nwave, GM_HLINE, 
				 Memr[tfwhm])

	    case IS_HIST:
		call histplot (gp, Memr[wave], Memr[flux], nwave)

		if (and (tbar, HAS_VERT) > 0)
		    call pltbar (gp, Memr[wave], Memr[flux], nwave, GM_VLINE,
				 Memr[terror])
	    }

	    flux = flux + nwave
	    err = err + nwave
	    fwhm = fwhm + nwave
	}

end
