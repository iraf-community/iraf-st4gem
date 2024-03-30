# COUNTRATE -- Compute the count rate of the HST

#* HISTORY *
#* B.Simon	14-Apr-93	original
#* B.Simon	25-Oct-93	modified because wfp2 != wfpc2
#* B.Simon	07-Jan-94	strdic() replaced by word_match()
#* B.Simon	15-Apr-94	resample on wavelength set after computations
#* B.Simon	31-May-94	rewritten for new version of synphot
#* B.Simon	02-Feb-95	handle case where instrument is a list
#* B.Simon	27-Jun-96	renormalize to vega spectrum
#* B.Simon	18-Jul-96	convert to output units with setform
#* B.Simon	06-Aug-96	use spectrum as wavelength set, if file
#* B.Simon	14-Aug-96	added magform as task parameter
#* B.Simon	21-Nov-96	added redlaw as task parameter
#* B.Simon	16-Oct-97	added acs to instrument list
#* B.Simon	13-Oct-98	resample spectrum for all instruments
#* B.Simon	21-Sep-00	compute pivlam & fwhmlam before resampling

procedure countrate

#--
pointer	spectrum	# spectrum to calculate
pointer	magnitude	# magnitude and passband of spectrum
pointer	instrument	# science instrument
pointer	detector	# detector used
pointer	spec_elem	# spectral elements used
pointer	aperture	# aperture / field of view
real	cenwave		# central wavelength
real	exptime		# exposure time in seconds
real	reddening	# interstellar reddening E(B-V)
pointer	redlaw		# reddening law used (gal1-3, smc, lmc, xgal)
pointer	output		# name of output table
pointer	form		# form (units) of output
pointer	magform		# form (units) of magnitude
pointer	wavecat		# catalog of wavelength tables
real	refwave		# wavelength at which counts are measured
real	flux_tot	# estimated total flux (output)
real	flux_ref	# estimated flux at reference wavelength (output)
bool	verbose		# print results to STDOUT ?
real	hstarea		# telescope area in cm^2
pointer	graphtab	# graph table
pointer	comptab		# component lookup table

int	nwave, degree
pointer	sp, key, mode, wavetab, newspec, pcode, ocode, scode
pointer	wave, band, flux
real	pivot, fwhm

string	blank        ""
string	notband      "Not a passband"
string	notspec      "Not a spectrum"

bool	clgetb(), rdwavecat()
int	tbtacc()
real	clgetr(), pivlam(), fwhmlam(), effstim(), getcounts()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (spectrum, SZ_FNAME, TY_CHAR)
	call salloc (magnitude, SZ_FNAME, TY_CHAR)
	call salloc (instrument, SZ_FNAME, TY_CHAR)
	call salloc (detector, SZ_FNAME, TY_CHAR)
	call salloc (spec_elem, SZ_FNAME, TY_CHAR)
	call salloc (aperture, SZ_FNAME, TY_CHAR)
	call salloc (redlaw, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (form, SZ_FNAME, TY_CHAR)
	call salloc (magform, SZ_FNAME, TY_CHAR)
	call salloc (wavecat, SZ_FNAME, TY_CHAR)
	call salloc (graphtab, SZ_FNAME, TY_CHAR)
	call salloc (comptab, SZ_FNAME, TY_CHAR)

	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (mode, SZ_FNAME, TY_CHAR)
	call salloc (wavetab, SZ_FNAME, TY_CHAR) 
	call salloc (newspec, SZ_LINE, TY_CHAR)
	call salloc (pcode, 2*SZ_LINE, TY_INT)

	ocode = pcode
	scode = pcode + SZ_LINE

	# Read task parameters

	call clgstr ("spectrum", Memc[spectrum], SZ_FNAME)
	call clgnone ("magnitude", Memc[magnitude], SZ_FNAME)
	call clgstr ("instrument", Memc[instrument], SZ_FNAME)

	call clgstr ("detector", Memc[detector], SZ_FNAME)
	call clgstr ("spec_elem", Memc[spec_elem], SZ_FNAME)
	call clgstr ("aperture", Memc[aperture], SZ_FNAME)

	cenwave = clgetr ("cenwave")
	exptime = clgetr ("exptime")
	reddening = clgetr ("reddening")
	call clgstr ("redlaw", Memc[redlaw], SZ_FNAME)

	call clgnone ("output", Memc[output], SZ_FNAME)
	call clgstr ("form", Memc[form], SZ_FNAME)
	call clgstr ("magform", Memc[magform], SZ_FNAME)
	call clgstr ("wavecat", Memc[wavecat], SZ_FNAME)

	refwave = clgetr ("refwave")
	verbose = clgetb ("verbose")

	hstarea = clgetr ("area")
	call clgstr ("grtbl", Memc[graphtab], SZ_FNAME)
	call clgstr ("cmptbl", Memc[comptab], SZ_FNAME)

	call put_hstarea (hstarea)

	# Check the form parameters

	call chkform (Memc[form])
	call chkform (Memc[magform])

	# Build and compile the commands to compute throughput and spectrum

	call inisyntab
	call buildmode (Memc[instrument], Memc[detector], Memc[spec_elem],
			Memc[aperture], Memc[mode], SZ_LINE)

	call syncompile (Memc[mode], Memi[ocode], SZ_LINE)

	call buildspec (Memc[spectrum], Memc[magnitude], Memc[magform],
			exptime, reddening, Memc[redlaw], Memc[newspec], 
			SZ_LINE)

	call syncompile (Memc[newspec], Memi[scode], SZ_LINE)

	# Compute the wavelength set

	if (tbtacc (Memc[spectrum]) == NO) {
	    call getwavelen (blank, Memc[graphtab], Memc[comptab], 
			     Memi[pcode], 2, SZ_LINE, wave, nwave)
	} else {
	    call getwavelen (Memc[spectrum], Memc[graphtab], Memc[comptab], 
			     Memi[pcode], 2, SZ_LINE, wave, nwave)
	}

	call malloc (band, nwave, TY_REAL)
	call malloc (flux, nwave, TY_REAL)

	# Compute the throughput and spectrum

	call syncalc (Memi[ocode], SZ_LINE, NULL, nwave, Memr[wave], 
		      Memc[graphtab], Memc[comptab], Memr[band], degree)

	if (degree != 0)
	    call printerr_str (notband, Memr[band])

	call syncalc (Memi[scode], SZ_LINE, NULL, nwave, Memr[wave], 
		      Memc[graphtab], Memc[comptab], Memr[flux], degree)

	if (degree != 1)
	    call printerr_str (notspec, Memr[flux])

	# Compute bandpass statistics

	pivot = pivlam (nwave, Memr[wave], Memr[band])
	fwhm = fwhmlam (nwave, Memr[wave], Memr[band])

	# Resample bandpass and flux on a wavelength set appropriate
	# for the instrument and grating combination

	Memc[wavetab] = EOS

	if (rdwavecat (Memc[wavecat], Memc[mode], Memc[wavetab], SZ_FNAME))
	    call resamp (Memc[wavetab], Memc[mode], cenwave, 
			 nwave, wave, band, flux)

	# Compute total counts and counts at specified wavelength

	flux_tot = effstim (nwave, Memr[wave], Memr[band], 
			     Memr[flux], Memc[form])

	call amulr (Memr[band], Memr[flux], Memr[flux], nwave)
	call setoutform  (Memc[form], Memc[wavetab], nwave, 
			  Memr[wave], Memr[flux])

	flux_ref = getcounts (Memc[form], refwave, nwave, 
			       Memr[wave], Memr[flux])

	# Write results to parameter file and STDOUT

	call clputr ("flux_tot", flux_tot)
	call clputr ("flux_ref", flux_ref)

	if (verbose)
	    call wrt_counts (Memc[mode], Memc[newspec], Memc[form],
			     refwave, pivot, fwhm, flux_tot, flux_ref)

	call putcounts (Memc[output], Memc[form], Memc[graphtab],
			Memc[comptab], Memc[mode], Memc[newspec], 
			exptime, nwave, Memr[wave], Memr[flux])

	# Free allocated memory

	call clssyntab

	call mfree (wave, TY_REAL)
	call mfree (band, TY_REAL)
	call mfree (flux, TY_REAL)
	call sfree (sp)

end
