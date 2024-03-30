#* HISTORY *
#* B.Simon	27-May-94	Original

# BUILDMODE -- Build an observation mode string from the individual components

procedure buildmode (instrument, detector, spec_elem, aperture, mode, maxch)

char	instrument[ARB]	# i: science instrument
char	detector[ARB]	# i: detector used
char	spec_elem[ARB]	# i: spectral elements used
char	aperture[ARB]	# i: aperture / field of view
int	mode[ARB]	# o: o
int	maxch		# i: maximum length of pseudocode
#--
pointer	sp, temp

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (temp, maxch, TY_CHAR)

	# Build mode string from component keywords

	call sprintf (Memc[temp], maxch, "%s %s %s %s")
	call pargstr (instrument)
	call pargstr (detector)
	call pargstr (spec_elem)
	call pargstr (aperture)

	call strlwr (Memc[temp])

	# Convert to standard form 

	call exp_rewrite (Memc[temp], mode, maxch)

	call sfree (sp)

end
