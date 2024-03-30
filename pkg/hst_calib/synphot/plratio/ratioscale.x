include	<tbset.h>
include "../plspec/plspec.h"

# RATIOSCALE -- Scale the ratio from plratio
#
# Sep 1993 H.A.Bushouse: Resample sphot data to same waveset as spec data
# Nov 1993 Bernie Simon; build ebmv string with applyebmv()

procedure ratioscale (nwave, wave, speclist, nspec, sphotlist, nsphot, 
	              phratio, sigma, nphot, ebmv1, ebmv2, nebmv, form,
	              xmin, xmax, ymin, ymax,
	              rms, bias, chisqr, nsum) 

int	nwave			# i: number of wavelenths
pointer	wave			# i: pointer to wavelength array
char	speclist[SZ_FNAME,ARB]	# i: list of synthetic spectrum files
int	nspec			# i: number of spectra
char	sphotlist[SZ_FNAME,ARB] # i: list of spectrophotometry files
int	nsphot			# i: number of spectrophotometry spectra
real	phratio[ARB]		# i: array of photometry/synphot points
real	sigma[ARB]		# i: 1-sigma error bars for phratio
int	nphot			# i: number of photometry ratios
real	ebmv1			# i: first ebmv value
real	ebmv2			# i: last ebmv value
int	nebmv			# i: number of ebmv values
char	form[ARB]		# i: form of data
real	xmin, xmax, ymin, ymax	# o: plot limits
real	rms, bias, chisqr, nsum	# o: statistical parameters from stats
#--
bool	status, mag
int	isphot, ispec, iebmv, iw, ic, nspwave
pointer sp, spec, spdat, spsig, spfwhm, spform, spmode, spstar, cmd
pointer sptwv, sptdt, sptsg
pointer	inform, ratio, grtbl, cmptbl, tot
real	debmv, ebmv, temp, xxmin, xxmax, yymin, yymax

int	strsearch()

begin
	# Allocate memory

	call smark (sp) 
	call salloc (spec, nwave, TY_REAL)
	call salloc (ratio, nwave, TY_REAL)
	call salloc (spdat, nwave, TY_REAL)
	call salloc (spsig, nwave, TY_REAL)
	call salloc (spform, SZ_COLUNITS, TY_CHAR)
	call salloc (spmode, SZ_LINE, TY_CHAR)
	call salloc (spstar, SZ_FNAME, TY_CHAR)
	call salloc (sptwv, MAXWAVE, TY_REAL)
	call salloc (sptdt, MAXWAVE, TY_REAL)
	call salloc (sptsg, MAXWAVE, TY_REAL)
	call salloc (spfwhm, MAXWAVE, TY_REAL)
	call salloc (cmd, SZ_LINE, TY_CHAR)
	call salloc (inform, SZ_COLUNITS, TY_CHAR)
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR)

	# Check for magnitudes in form string

	mag = strsearch(form,"mag") > 0 || strsearch(form,"MAG") > 0

	debmv = (ebmv2 - ebmv1) / max (1,nebmv-1)
	# Save the user limits

	xxmin = xmin
	xxmax = xmax
	yymin = ymin
	yymax = ymax

	# Initialize limits

	if (IS_INDEFR(xmin) || IS_INDEFR(xmax))
	    call xlimit (Memr[wave], nwave, xmin, xmax)

	if (mag) {
	    ymin = 0.0
	    ymax = 0.0
	} else {
	    ymin = 1.0
	    ymax = 1.0
	}

	# Get graph and component table names

	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)


	do isphot = 1, nsphot {
	    nspwave = MAXWAVE
	    call loadsphot (nspwave, sphotlist[1,isphot], Memr[sptwv], 
	                   Memr[sptdt], Memr[sptsg], Memr[spfwhm], 
	                   Memc[spform], Memc[spmode], Memc[spstar]) 

	    # Resample the spectrum and sigmas
	    call linterp (Memc[spform], nspwave, Memr[sptwv], Memr[sptdt],
			  nwave, Memr[wave], Memr[spdat])

	    call linterp (Memc[spform], nspwave, Memr[sptwv], Memr[sptsg],
			  nwave, Memr[wave], Memr[spsig])

	    call sphotform (nwave, Memr[wave], Memr[spdat], Memr[spsig], 
			    Memc[spform], Memr[spdat], Memr[spsig], form) 

	    # Loop over spectra

	    do ispec = 1, nspec {

		# Loop over E(B-V) values

		do iebmv = 1, max (1, nebmv) {
		    ebmv = ebmv1 + (iebmv-1)*debmv
		    call applyebmv (speclist[1,ispec], Memc[cmd], 
				    ebmv, SZ_LINE)

		    iw = 1
		    call compspec (Memc[cmd], iw, Memc[grtbl], Memc[cmptbl],
				   nwave, wave, spec, Memc[inform]) 
		    call specform (nwave, Memr[wave], Memr[spec], Memc[inform],
				   Memr[spec], form, status) 

		    do ic = 1, nwave {
			if (Memr[spec+ic-1] == 0 || 
			    Memr[spec+ic-1] == INDEFR ||
			    Memr[spdat+ic-1] == INDEFR) {

			    Memr[ratio+ic-1] = INDEFR

			} else if ( mag ) {
			    Memr[ratio+ic-1] = 
					    Memr[spdat+ic-1] - Memr[spec+ic-1]

			} else {
			    Memr[ratio+ic-1] = 
					    Memr[spdat+ic-1] / Memr[spec+ic-1]
			}
		    }

		    call xylimit (Memr[wave], Memr[ratio], nwave, 
				  xmin, xmax, ymin, ymax)
		}
	    }
	}

	# Check the photometry ratio

	if (nphot > 0) {

	    # Allocate memory for sum and difference

	    call salloc( tot, nphot, TY_REAL )

	    # Add 1-sigma errors to data before taking limit

	    do ic = 1, nphot {
		if (! IS_INDEFR (sigma[ic]))
		    Memr[tot+ic-1] = phratio[ic] + sigma[ic]
		else
		    Memr[tot+ic-1] = phratio[ic]
	    }
	    call ylimit (Memr[tot], nphot, ymin, ymax) 

	    # Subtract 1-sigma errors from data before taking limit

	    do ic = 1, nphot {
		if (! IS_INDEFR (sigma[ic]))
		    Memr[tot+ic-1] = phratio[ic] - sigma[ic]
		else
		    Memr[tot+ic-1] = phratio[ic]
	    }
	    call ylimit (Memr[tot], nphot, ymin, ymax) 
	}

	# Reset the plot limits to the user limits where they are not INDEF
	if (! IS_INDEFR (xxmin))
	    xmin = xxmin

	if (! IS_INDEFR (xxmax))
	    xmax = xxmax

	if (! IS_INDEFR (yymin))
	    ymin = yymin

	if (! IS_INDEFR (yymax))
	    ymax = yymax

	if (mag) {
	    temp = ymin
	    ymin = ymax
	    ymax = temp
	}

	# Get statistics for ratio

	call stats (Memr[spdat], Memr[spsig], Memr[spec], nwave, form,
	               chisqr, bias, rms, nsum)

	call sfree (sp)
end
