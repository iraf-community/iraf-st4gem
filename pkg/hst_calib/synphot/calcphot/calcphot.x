#* HISTORY *
#* B.Simon	17-May-94	Original
#* B.Simon	20-Sep-00	Updated calculation of wavelength set
#* B.Simon	30-Oct-00	Added call to exp_rewrite before syncompile

# CALCPHOT -- Calculate synthetic photometry for a set of obsmodes and spectra

procedure calcphot ()

#--
pointer	obsmode		# File containing list of modes
pointer	spectrum	# File containing list of spectra
pointer form		# Form of output synthetic data
pointer func		# Function used to calculate output
pointer	vzero		# List of values for variable zero
pointer	wavetab		# File containing wavelength array
pointer output		# Optional output table name
bool	append		# Append results to existing output table?
pointer	grtbl		# Instrument graph table
pointer cmptbl		# Instrument component table
real	hstarea		# Telescope area
real	result		# Result of synthetic photometry

bool	first
extern	getsynvar
int	nwave, nband, irow, ispec, degree
pointer	sp, mode, mode1, mode2, spec, pcode, ocode, ocode1, ocode2, scode
pointer	wave, band1, band2, flux, optr, sptr
real	v0

string	nomode   "No obsmode. Must set form to counts or obmag"
string	notspec  "Not a spectrum"

bool	clgetb()
int	is_count()
pointer	rdlist(), nxtlist(), nxtvzero(), locpr()
real	clgetr()

begin

	# Allocate memory for strings and arrays

	call smark (sp)
	call salloc (obsmode, SZ_FNAME, TY_CHAR)
	call salloc (spectrum, SZ_FNAME, TY_CHAR)
	call salloc (form, SZ_FNAME, TY_CHAR)
	call salloc (func, SZ_FNAME, TY_CHAR)
	call salloc (vzero, SZ_FNAME, TY_CHAR)
	call salloc (wavetab, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR )
	call salloc (grtbl, SZ_FNAME, TY_CHAR)
	call salloc (cmptbl, SZ_FNAME, TY_CHAR) 

	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (mode1, SZ_LINE, TY_CHAR)
	call salloc (mode2, SZ_LINE, TY_CHAR)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (pcode, 4*SZ_LINE, TY_INT)

	scode = pcode
	ocode = pcode + SZ_LINE
	ocode1 = pcode + 2 * SZ_LINE
	ocode2 = pcode + 3 * SZ_LINE

	# Read in task parameters 

	call clgnone ("obsmode", Memc[obsmode], SZ_FNAME)
	call clgnone ("spectrum", Memc[spectrum], SZ_FNAME) 
	call clgstr ("form", Memc[form], SZ_FNAME) 
	call clgstr ("func", Memc[func], SZ_FNAME) 
	call clgstr ("vzero", Memc[vzero], SZ_FNAME)
	call clgnone ("wavetab", Memc[wavetab], SZ_FNAME)
	call clgnone ("output", Memc[output], SZ_FNAME)
	append = clgetb ("append") 

	call clgstr ("grtbl", Memc[grtbl], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptbl], SZ_FNAME)
	hstarea = clgetr ("area")

	call put_hstarea (hstarea)

	# Check form parameter

	call chkform (Memc[form])

	# Create output table

	call mkphottab (Memc[output], append, Memc[func], irow)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Check for missing obsmode, substitute one

	if (Memc[obsmode] == EOS) {
	    if (is_count (Memc[form]) == NO)
		call printerr_str (nomode, Memc[form])

	    call strcpy ("1.0", Memc[obsmode], SZ_LINE)
	}

	# Process the obsmode, spectrum,  and vzero strings

	optr = rdlist (Memc[obsmode])
	sptr = rdlist (Memc[spectrum])
	call rdvzero (Memc[vzero])

	# Process each observation mode / bandpass

	first = true
	result = 0.0

	while (nxtlist (optr, Memc[mode], SZ_LINE) != EOF) {

	    # Compile the observation mode expression

	    call exp_rewrite (Memc[mode], Memc[mode1], SZ_LINE)
	    call syncompile (Memc[mode1], Memi[ocode], SZ_LINE)

	    # Delay computing results until wavelength set is determined
	    # This can only be done after the first spectrum is read

	    if (! first) {
		call photband (Memc[grtbl], Memc[cmptbl], Memc[mode],
			       nwave, Memr[wave], SZ_LINE, nband, 
			       Memc[mode1], Memc[mode2], Memi[ocode1], 
			       Memi[ocode2], Memr[band1], Memr[band2])
	    }

	    # Process each target / spectrum

	    while (nxtlist (sptr, Memc[spec], SZ_LINE) != EOF) {
		call syncompile (Memc[spec], Memi[scode], SZ_LINE)

		# The first obsmode is processed here since we have to wait
		# until first spectrum is read to compute the wavelength set

		if (first) {
		    first = false

		    # Compute wavelength set and allocate arrays

		    call getwavelen (Memc[wavetab], Memc[grtbl], Memc[cmptbl], 
				     Memi[pcode], 2, SZ_LINE, wave, nwave)

		    call salloc (band1, nwave, TY_REAL)
		    call salloc (band2, nwave, TY_REAL)
		    call salloc (flux, nwave, TY_REAL)

		    call photband (Memc[grtbl], Memc[cmptbl], Memc[mode],
				   nwave, Memr[wave], SZ_LINE, nband, 
				   Memc[mode1], Memc[mode2], Memi[ocode1], 
				   Memi[ocode2], Memr[band1], Memr[band2])
		}

		# Process each value of vzero

		ispec = 0
		while (nxtvzero (v0) != EOF) {
		    if (IS_INDEFR(v0))
			v0 = 0.0

		    # Compute flux of target /spectrum

		    call syncalc (Memi[scode], SZ_LINE, locpr (getsynvar),
				  nwave, Memr[wave], Memc[grtbl], 
				  Memc[cmptbl], Memr[flux], degree)

		    if (degree != 1)
			call printerr_str (notspec, Memc[spec])

		    # Compute requested function of bandpass and spectrum

		    call photstim (Memc[form], Memc[func], nwave, Memr[wave], 
				   nband, Memr[band1], Memr[band2], 
				   Memr[flux], result)

		    # Print results to STDOUT and save in table

		    call photoutput (Memc[form], Memc[func], Memc[mode1], 
				     Memc[mode2],  Memc[spec], nband, v0, 
				     result, ispec)

		    call wrtphottab (Memc[form], Memc[mode], Memc[spec],
				     result, irow)
				      
		}
	    }
	}

	# Write result of last calculation to result

	call clputr ("result", result)

	# Free structures and memory

	call clsphottab
	call clssyntab

	call freelist (optr)
	call freelist (sptr)

	if (! first)
	    call mfree (wave, TY_REAL)

	call sfree (sp)
end

