# SETNULL -- Set outlier values to INDEF

#* HISTORY *
#* B.Simon	06-Apr-94	first code

procedure setnull (form, nwave, wave, spec)

char	form[ARB]	# i: flux form name
int	nwave		# i: number of array elements
real	wave[ARB]	# i: wavelength set, in angstroms
real	spec[ARB]	# u: flux array
#--
int	iwave
int	is_magunit()

begin
	# In non-magnitude units zero stands as a proxy for INDEF
	# In magnitude units one hundred stands as a proxy
	# Restore the INDEFs for software that expects to find them

	if (is_magunit (form) == NO) {
	    do iwave = 1, nwave {
		if (spec[iwave] <= 0.0)
		    spec[iwave] = INDEFR
	    }

	} else {
	    do iwave = 1, nwave {
		if (spec[iwave] >= 100.0)
		    spec[iwave] = INDEFR
	    }
	}
end
