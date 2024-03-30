# SPHOTFORM -- Convert the form of spectro-photometry data and errors.

procedure sphotform( nwave, wave, datold, sigold, formold, 
	             datnew, signew, formnew )

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: Wavelength array
real	datold[ARB]	# i: Old data
real	sigold[ARB]	# i: Old 1 sigma error bars
char	formold[ARB]	# i: Form of input data
real	datnew[ARB]	# o: Output data
real	signew[ARB]	# o: Output sigmas
char	formnew[ARB]	# o: Form of output data

int	i, nscan
int	strsearch(), strlen(), scan()
bool	status
bool	streq()

# feb 1989 keith horne @ stsci - original version
# Oct 1989 Dave Bazell - SPP version

begin

	# return if no conversion is required

	if( nwave <= 0 )
	   return

	if ( strlen(formnew) <= 0 ) {
	   call printf("Output form for SPHOTFORM? ")
	   nscan = scan()
	      call gargstr(formnew)
	}

# maybe no conversion is necessary

	if (streq(formold, formnew)) {
	   do i=1,nwave {
	      datnew[i] = datold[i]
	      signew[i] = sigold[i]
	   }
	   return
	}

# make new array of data and 1-sigma brighter than data

	if( strsearch(formold, "mag" ) > 0 ) {

	   do i=1,nwave {

	      datnew[i] = datold[i]

	      if( !IS_INDEFR (datold[i]) && !IS_INDEFR (sigold[i]) ) {
	         if( sigold[i] > 0. )
	            signew[i] = datold[i] - sigold[i]	# brighter by 1 sigma
	         else
	            signew[i] = datold[i] + sigold[i]	# fainter by 1 sigma

	      } else
	         signew[i] = INDEFR
	   }

	} else {

	   do i=1,nwave {

	      datnew[i] = datold[i]

	      if( !IS_INDEFR (datold[i]) && !IS_INDEFR (sigold[i]) ) {
	         if( sigold[i] > 0. )
	            signew[i] = datold[i] + sigold[i]	# brighter by 1 sigma
	         else
	            signew[i] = datold[i] - sigold[i]	# fainter by 1 sigma

	      } else
	         signew[i] = INDEFR
	   }
	}

# transform both new arrays

	call specform( nwave, wave, datnew, formold, datnew, formnew,status )
	call specform( nwave, wave, signew, formold, signew, formnew,status )

# take difference to get 1-sigma of transformed data

	if( strsearch( formnew, "mag" ) > 0 ) {
	   do i=1,nwave {
	      if ( !IS_INDEFR (datnew[i]) && !IS_INDEFR (signew[i]))
	         signew[i] = datnew[i] - signew[i]
	   }

	} else {
	   do i=1,nwave {
	      if ( !IS_INDEFR (datnew[i]) && !IS_INDEFR (signew[i]) )
	         signew[i] = signew[i] - datnew[i]
	   }
	}
end
