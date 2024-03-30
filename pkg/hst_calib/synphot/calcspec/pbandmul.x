# PBANDMUL -- Multiply a spectrum by a passband

procedure pbandmul( script, grtbl, cmptbl, nwv, wv, form, iw, spec )

char	script[ARB]	# i: command string
char 	grtbl[ARB]	# i: Graph table name
char 	cmptbl[ARB]	# i: Component table name
int 	nwv		# i: number of wavelengths
pointer	wv		# i: pointer to wavelength set
char 	form[ARB]	# i: Form of spectrum, i.e. units
int	iw		# io: position in script
real 	spec[ARB]	# io: spectrum or spec * passband

int	i
int	strsearch()
pointer filt

begin

	# Evaluate instrument passband
	call compband( script, iw, grtbl, cmptbl, nwv, wv, filt)

	# If spectrum is in magnitudes convert filt to mags and add
        if( strsearch(form,"mag") > 0 )
           do i=1,nwv {
              if( Memr[filt+i-1] > 0. && Memr[filt+i-1] != INDEFR &&
	                                 !IS_INDEFR (spec[i]) ) 
	         spec[i] = spec[i] - 2.5 * alog10( Memr[filt+i-1] )
              else
	         spec[i] = INDEFR
	   }

 	# If spectrum not magnitudes, just multiply
        else
	   call bmulr ( spec, Memr[filt], spec, nwv )

	call mfree(filt, TY_REAL)

end
