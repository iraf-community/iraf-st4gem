# FILLNULL -- Replace nulls with interpolated valaues

int procedure fillnull (fillval, nwave, wave, data)

real	fillval		# i: value used to replace nulls at ends
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	data[ARB]	# u: data array
#--
int	inil, iwave, jwave, status
real	frac

begin
	# inil is the index of the start of the most recent run of nulls

	inil = 0
	do iwave = 1, nwave {
	    if (inil == 0) {

		# Mark the start of a run of null values
		if (IS_INDEFR(data[iwave]))
		    inil = iwave

	    } else {
		if (! IS_INDEFR(data[iwave])) {

		    # Fill in the run when a run of nulls is found
		    if (inil == 1) {

			# Set data to fillval if there is no bracket
			do jwave = inil, iwave - 1
			    data[jwave] = fillval

		    } else {

			# Do linear interpolation if there is a bracket
			do jwave = inil, iwave - 1 {
			    frac = (wave[jwave] - 
				    wave[inil-1])/ (wave[iwave] - wave[inil-1])
			    data[jwave] = data[inil-1] * (1.0 - frac) + 
					  data[iwave] * frac
			}
		    }
		    inil = 0
		}
	    }
	}

	# Final run of nulls has no bracket, so set to fillval
	if (inil != 0) {
	    do jwave = inil, nwave
		data[jwave] = fillval
	}

	# Return error status if all data values are null

	if (inil == 1) {
	    status = ERR
	} else {
	    status = OK
	}
	return (status)
end
