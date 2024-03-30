# PHOTFORM -- Convert the form of photometric data and errors.

procedure photform( mode, ndata, nwave, wave, datold, sigold, formold,
	                                      datnew, signew, formnew )

# input:
char	mode[ARB]	# i: passband
int	ndata		# i: number of data points
int	nwave		# i: number of wavelength points
pointer	wave		# i: pointer to wavelength array
real	datold[ARB]	# i: old data
real	sigold[ARB]	# i: old 1-sigma error bars
char	formold[ARB]	# i: form of old data on input (e.g. 'fnu')
char	formnew[ARB]	# i: form of spectrum desired on output
# output:
real	datnew[ARB]	# o: output data transformed to new form
real	signew[ARB]	# o: output sigmas transformed to new form

# jun 1990 keith horne @ stsci - original version
# Oct 1990 Dave Bazell  Convert to spp
# Nov 1993 Bernie Simon Modify for new expression evaluator

bool	magold, magnew
bool	streq()
int	ic
int	strsearch()
real	band[1], area, rateold, ratenew, scale, sigma, data, bright, add
char	cmd[SZ_LINE], grtbl[SZ_FNAME], cmptbl[SZ_FNAME]

begin

	# Return if no conversion is required
	if ( ndata <= 0 ) return

	# Maybe no conversion is necessary

	if( streq(formold,formnew) ) {
	   do ic = 1,ndata {
	      datnew[ic] = datold[ic]
	      signew[ic] = sigold[ic]
	   }
	   return
	}

	# Get the graph table and component table names

	call clgstr( "grtbl", grtbl, SZ_FNAME)
	call clgstr( "cmptbl", cmptbl, SZ_FNAME)

	call strlwr( formold )
	call strlwr( formnew )

	magold = false
	magnew = false
	if (strsearch( formold, "mag" ) > 0 )
	   magold = true
	if (strsearch( formnew, "mag" ) > 0 )
	   magnew = true

	# count rate for old units
	if( strsearch( formold,"counts") > 0 ||
	    strsearch( formold,"obmag") > 0 ) {
	   call get_hstarea (area)
	   rateold = 1.0 / area

	# construct string "unit(0,form)" (e.g. unit(0,abmag) )
	} else if( magold ) {
	   call sprintf (cmd, SZ_LINE, "unit(0,%s)")
	   call pargstr (formold)
	   call synphot( nwave, wave, band, cmd, "counts", mode, 
	                 grtbl, cmptbl, rateold )

	# construct string "unit(1,form)" (e.g. unit(1,flam) )
	} else {
	   call sprintf (cmd, SZ_LINE, "unit(1,%s)")
	   call pargstr (formold)
	   call synphot( nwave, wave, band, cmd, "counts", mode, 
	                 grtbl, cmptbl, rateold )
	}

	# Count rate for new units
	if( strsearch( formnew,"counts") > 0 ||
	    strsearch( formnew,"obmag") > 0 )
	   ratenew = 1.0 / area

	# construct string "unit(0,form)" (e.g. unit(0,abmag) )
	else if( magnew ) {
	   call sprintf (cmd, SZ_LINE, "unit(0,%s)")
	   call pargstr (formnew)
	   call synphot( nwave, wave, band, cmd, "counts", mode, 
	                 grtbl, cmptbl, ratenew )

	# construct string "unit(1,form)" (e.g. unit(1,flam) )
	} else {
	   call sprintf (cmd, SZ_LINE, "unit(1,%s)")
	   call pargstr (formnew)
	   call synphot( nwave, wave, band, cmd, "counts", mode, 
	                 grtbl, cmptbl, ratenew )
	}

	scale = rateold / ratenew

	# loop thru data points
	do ic=1,ndata {

	   # Data and 1-sigma brighter than data
	   data = datold[ic]
	   sigma = sigold[ic]
	   if ( !IS_INDEFR (data) && !IS_INDEFR (sigma) ) {
	      if( magold )
	         bright = data - sigma
	      else 
	         bright = data + sigma
	   } else
	      bright = INDEFR

	   # transform
	   if( !magold && !magnew ) {
	      if ( !IS_INDEFR (data) )
	         data = data * scale
	      if (!IS_INDEFR (bright) )
	         bright = bright * scale

	   } else if( magold && magnew ) {
	      add = - 2.5 * alog10( scale )
	      if ( !IS_INDEFR (data) )
	         data = data + add
	      if (!IS_INDEFR (bright))
	         bright = bright + add

	   } else if( magold && !magnew ) {
	      add = - 2.5 * alog10( scale )
	      if (!IS_INDEFR (data) )
	         data = 10.** ( -0.4 * ( data + add ) )
	      if (!IS_INDEFR (bright))
	         bright = 10.** ( -0.4 * ( bright + add ) )

	   } else if( !magold && magnew ) {
	      add = - 2.5 * alog10( scale )
	      if (!IS_INDEFR (data) )
	         data = -2.5 * alog10( data ) + add
	      if (!IS_INDEFR (bright) )
	         bright = -2.5 * alog10( bright ) + add
	   }

	   # take difference to recover 1-sigma of transformed data

	   datnew[ic] = data
	   if ( magnew ) {
	      if ( !IS_INDEFR (data) && !IS_INDEFR (sigold[ic]) && !IS_INDEFR (bright))
	         signew[ic] = data - bright
	      else
	         signew[ic] = INDEFR
	   } else {
	      if ( !IS_INDEFR (data) && !IS_INDEFR (sigold[ic]) && !IS_INDEFR (bright))
	         signew[ic] = bright - data
	      else
	         signew[ic] = INDEFR
	   }

	}

end
