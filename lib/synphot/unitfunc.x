# UNITFUNC -- Comput a constant spectrum with specified units

procedure unitfunc (const, units, nwave, wave, spec)

real	const		# i: Constant value
char	units[ARB]	# i: Physical units of constant
int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: output spectrum, in photlam
#--
int	done

string	badunits  "Unknown flux units for function unit"

int	anytophot()
errchk	synphoterr

begin
	# Copy to spectrum and convert to photlam

	call amovkr (const, spec, nwave)

	iferr {
	    done = anytophot (units, nwave, wave, spec)
	} then {
	    done = NO
	} else {
	    done = YES
	}

	if (done == NO)
	    call synphoterr (badunits, units)

end
