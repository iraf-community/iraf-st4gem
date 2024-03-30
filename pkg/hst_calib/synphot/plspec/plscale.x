include	"plspec.h"

# PLSCALE -- Autoscale plot given input spectrum files

# Oct 1989 Dave Bazell; Modularized and converted to SPP
# Dec 1989 DB; Fix bug in E(B-V) loop
# Jan 1990 DB; Insert 'ebmv' command before '*' command
# Oct 1990 DB; use photform to chage form of photometry
# Nov 1990 DB; Fix scaling of photometry data
# Aug 1993 Howard Bushouse; Fixed bug in scaling spectral data in obmag units
# Sep 1993 HB; Add warning for 2 mode photometry data
# Nov 1993 Bernie Simon; build ebmv string with applyebmv()

procedure plscale(modelist, speclist, splist, plist, form, nwave, wave,
	          grtbl, cmptbl, ebmv1, ebmv2, nebmv, hstarea, nmode, nspec,
	          nsphot, nphot, xmin, xmax, ymin, ymax)

char	modelist[SZ_LINE,ARB]	# i: Passband list
char	speclist[SZ_LINE,ARB]	# i: Spectrum list
char	splist[SZ_LINE,ARB]	# i: Spectrophotometry list
char	plist[SZ_LINE,ARB]	# i: Photometry list
char	form[ARB]	# i: Spectrum form
int	nwave		# i: Number of wavelengths in wave
pointer	wave		# i: Pointer to wavelength array
char	grtbl[ARB]	# i: Instrument graph table
char	cmptbl[ARB]	# i: Instruemnt component table
real	ebmv1		# i: First E(B-V) value
real	ebmv2		# i: Last E(B-V) value
int	nebmv		# i: Number of E(B-V) values
real	hstarea		# i: HST area
int	nmode		# o: Number of modes
int	nspec		# o: Number of spectra
int	nsphot		# o: Number of spectrophotometry files
int	nphot		# o: Number of photometry files
real	xmin		# o: Minimum x value of plot
real	xmax		# o: Maximum x value of plot
real	ymin		# o: Minimum y value of plot
real	ymax		# o: Maximum y value of plot
#--
bool	status, mag
char	mode[SZ_LINE], mode1[SZ_LINE], mode2[SZ_LINE]
char	pform[SZ_FNAME,MAXPHOT], pmode[SZ_LINE,MAXPHOT]
int	idat, npdat, nband, ip, isp, nspwave, imode
int	iebmv, ispec, ib, iw
real	piv, fwhm, dat, sig, ebmv, temp
real	xxmin, xxmax, yymin, yymax, ydif
pointer	sp, pdat, psig, pstar, spwave, spdat, spsig, spfwhm, inform, spmode
pointer	spstar, specname, spec, filt, flterr, command, errmsg, band

string	novertical  "No vertical extent in plot:  %g - %g"
string	bandfmt	    "(%s)*band(%s)"
string	expfmt	    "(%s)*(%s)"
string	nulcmd	    ""

bool	is_simple()
int	strsearch()

begin
	
	# Allocate temporary memory
	
	call smark(sp)
	# Photometry
	call salloc (pdat, MAXPDAT, TY_REAL)
	call salloc (psig, MAXPDAT, TY_REAL)
	call salloc (pstar, SZ_FNAME, TY_CHAR)
	# Spectrophotometry
	call salloc (spwave, MAXWAVE, TY_REAL)
	call salloc (spdat, MAXWAVE, TY_REAL)
	call salloc (spsig, MAXWAVE, TY_REAL)
	call salloc (spfwhm, MAXWAVE, TY_REAL)	
	call salloc (inform, SZ_FNAME, TY_CHAR)
	call salloc (spmode, SZ_LINE, TY_CHAR)
	call salloc (spstar, SZ_FNAME, TY_CHAR)
	call salloc (specname, SZ_FNAME, TY_CHAR)
	call salloc (command, SZ_LINE, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)
	call salloc (band, MAXWAVE, TY_REAL)
	
	call printf("Autoscaling, %s.\n")
	call pargstr(form)
	
	# Check for magnitudes in form string
	
	call strlwr (form)
	mag = strsearch(form,"mag") > 0
	
	# Save the user limits
	
	xxmin = xmin
	xxmax = xmax
	yymin = ymin
	yymax = ymax
	
	# Initialize limits
	
	if (IS_INDEFR(xmin) || IS_INDEFR(xmax))
	    call xlimit (Memr[wave], nwave, xmin, xmax)

	ymin = INDEFR
	ymax = INDEFR

	# Autoscale the photometry data

	do ip=1,nphot {
	    npdat = MAXPDAT
	    call loadphot (npdat, plist[1,ip], Memr[pdat], Memr[psig], 
			   pform, pmode, Memc[pstar])

	    do idat=1,npdat {
		call splitmode (pmode[1,idat], nband, mode1, mode2)

		if (nband == 1)  {
		    call synphot (nwave, wave, Memr[band], nulcmd, "pivlam", 
				  mode1, grtbl, cmptbl, piv)
		    call synphot (nwave, wave, Memr[band], nulcmd, "fwhmlam",
				  mode1, grtbl, cmptbl, fwhm)
		    
		    # Check photometry points one by one

		    dat = Memr[pdat+idat-1]
		    sig = Memr[psig+idat-1]

		    call photform (mode1, 1, nwave, wave, dat, sig, 
				   pform[1,idat], dat, sig, form)

		    xmin = min (xmin, piv - fwhm)
		    xmax = max (xmax, piv + fwhm)
		    
		    # If data and errors ok then use to get limits
		    # Just use data to get limits if errors are indef
		    
		    if (!IS_INDEFR(dat)) {
			if (IS_INDEFR(sig)) {
			    if (IS_INDEFR(ymin) || IS_INDEFR(ymax)) {
				ymin = dat
				ymax = dat
			    } else {
				ymin = min (ymin, dat)
				ymax = max (ymax, dat)
			    }

			} else {
			    if (IS_INDEFR(ymin) || IS_INDEFR(ymax)) {
				ymin = ymin - 2.0 * abs(sig)
				ymax = ymax + 2.0 * abs(sig)
			    } else {
				ymin = min (ymin, dat - 2.0*abs(sig))
				ymax = max (ymax, dat + 2.0*abs(sig))
			    }
			}
		    }
		} else {
	      call printf ("Cannot compute scale for 2 mode photometry data.\n")
	      call flush (STDOUT)
		}
	    }
	}
	
	# Add margin to ymin and ymax
	
	if (! IS_INDEFR(ymin) && ! IS_INDEFR(ymax)) {
	    ydif = ymax - ymin
	    ymin = ymin - 0.05 * ydif
	    ymax = ymax + 0.05 * ydif
	}
	
	# Autoscale the spectrophotometry
	
	do isp=1,nsphot {
	    nspwave = MAXWAVE
	    call loadsphot (nspwave, splist[1,isp], Memr[spwave], Memr[spdat], 
			    Memr[spsig], Memr[spfwhm], Memc[inform], 
			    Memc[spmode], Memc[spstar])
	    
	    call sphotform (nspwave, Memr[spwave], Memr[spdat], Memr[spsig],
			    Memc[inform], Memr[spdat], Memr[spsig], form)

	    call xylimit (Memr[spwave], Memr[spdat], nspwave,
			  xmin, xmax, ymin, ymax)
	    
	}	# next spectro-photometry file
	
	# Loop through extreme values of E(B-V).
	
	ebmv = ebmv1
	do iebmv=1,nebmv {
	    
	    # Loop through the spectra
	    do ispec=1,nspec {
		
		# Autoscale detected counts  (spectrum * passband * area)
		
		if ((strsearch(form,"counts") > 0 || 
		     strsearch(form,"photpix") > 0 ||
		     strsearch(form,"obmag") > 0) && nmode > 0)  {
		    do imode=1,nmode {
			call strcpy (modelist[1,imode], mode, SZ_LINE)
			call splitmode (mode, nband, mode1, mode2)

			do ib=1,nband {
			    if (is_simple (mode1)) {
				call sprintf(Memc[command], SZ_LINE, bandfmt)
			    } else {
				call sprintf(Memc[command], SZ_LINE, expfmt)
			    }
			    call pargstr (speclist[1,ispec])
			    call pargstr (mode1)

			    call applyebmv (Memc[command], Memc[command], 
					    ebmv, SZ_LINE)

			    iw = 1
			    call compspec (Memc[command],iw, grtbl, cmptbl, 
					   nwave, wave, spec, Memc[inform])
			    call specform (nwave, Memr[wave], Memr[spec], 
					   Memc[inform], Memr[spec], form, 
					   status)

			    call xylimit (Memr[wave], Memr[spec], nwave,
					  xmin, xmax, ymin, ymax)
			    call strcpy (mode2, mode1, SZ_LINE)
			}
		    }
		    
		    # Autoscale the spectrum
		    
		} else {
		    call applyebmv (speclist[1,ispec], Memc[command], 
				    ebmv, SZ_LINE)

		    iw = 1
		    call compspec (Memc[command], iw, grtbl, cmptbl, nwave,
				   wave, spec, Memc[inform])
		    call specform (nwave, Memr[wave], Memr[spec], Memc[inform],
				   Memr[spec], form, status)

		    call xylimit (Memr[wave], Memr[spec], nwave,
				  xmin, xmax, ymin, ymax)
		}
	    }	# next spectrum
	    ebmv = ebmv2
	}	# next e(b-v)
	
	# Reset the plot limits to the user limits where they are not INDEF
	
	if (! IS_INDEFR (xxmin))
	    xmin = xxmin
	
	if (! IS_INDEFR (xxmax))
	    xmax = xxmax
	
	if (! IS_INDEFR (yymin))
	    ymin = yymin
	
	if (! IS_INDEFR (yymax))
	    ymax = yymax
	
	# Fix up magnitude limits.  Check to see if form is mag and make
	# sure that there is really spectral data (nspec > 0)
	
	if (mag) {
	    temp = ymin
	    ymin = ymax
	    ymax = temp
	}
	
	if (ymin == ymax)  {
	    call sprintf(Memc[errmsg], SZ_LINE, novertical)
	    call pargr(ymin)
	    call pargr(ymax)
	    call error(1, Memc[errmsg])
	}
	
	call mfree (filt, TY_REAL)
	call mfree (flterr, TY_REAL)
	call sfree (sp)
end
