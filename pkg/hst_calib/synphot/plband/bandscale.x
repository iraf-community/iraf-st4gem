include	<tbset.h>

# BANDSCALE -- Autoscale real and synthetic instrument passbands

procedure bandscale(nband, bandlist, nwave, wave, grtbl, cmptbl, ylog,
	            xmin, xmax, ymin, ymax )

int	nband			# i: Number of bands in bandlist
char	bandlist[SZ_LINE, ARB]	# i: Synthetic passband commands
int	nwave			# i: Number of wavelengths
pointer	wave			# i: Pointer to wavelength array
char	grtbl[ARB]		# i: Instrument graph table
char	cmptbl[ARB]		# i: Instrument component table
bool	ylog			# i: Log of y data 
real	xmin			# u: min x value
real	xmax			# u: max x value
real	ymin			# u: min y value
real	ymax			# u: max y value
#--
int	ib, iw, nmode
pointer	band, sp, mode1, mode2
real	xxmin, xxmax, yymin, yymax

begin

	if (nband == 0)
	    return

	# Allocate temporary memory
	call smark (sp) 
	call salloc (mode1, SZ_LINE, TY_CHAR) 
	call salloc (mode2, SZ_LINE, TY_CHAR) 

	# Save the user limits
	xxmin = xmin
	xxmax = xmax
	yymin = ymin
	yymax = ymax

	# Set the x plot limits if INDEF
	if (IS_INDEFR(xmin) || IS_INDEFR(xmax))
	    call xlimit (Memr[wave], nwave, xmin, xmax)

	# Find appropriate scaling for the synthetic passbands
	do ib = 1, nband {

	    # split the passband into two modes if possible.
	    # Do the first mode.
	    call splitmode (bandlist[1,ib], nmode, Memc[mode1], Memc[mode2]) 

	    iw = 1
	    call compband (Memc[mode1], iw, grtbl, cmptbl, nwave, wave,
			   band) 

	   # Replace 0 and negatives with INDEF if plotting log
	    if (ylog)
		call toindef (Memr[band], nwave) 

	    call xylimit (Memr[wave], Memr[band], nwave, 
			  xmin, xmax, ymin, ymax) 
	    call mfree (band, TY_REAL)

	    # If band was split into two mode then do the second 
	    if (nmode == 2) {
		iw = 1
		call compband (Memc[mode2], iw, grtbl, cmptbl, 
			       nwave, wave, band) 

		# Replace 0 and negatives with INDEF if plotting log
		if ( ylog )
		    call toindef (Memr[band], nwave) 

		call xylimit (Memr[wave], Memr[band], nwave, 
			      xmin, xmax, ymin, ymax) 
		call mfree (band, TY_REAL)
	    }
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

	call sfree (sp)
end
