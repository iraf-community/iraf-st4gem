define	MAXPHOT		13

#* HISTORY *
#* B.Simon	06-Oct-94	original
#* B.Simon	01-Mar-00	added wpeak

# COMPPAR -- Compute bandpass parameters and write to table

procedure comppar (obsmode, output, photlist, refwave,
		   wavtable, grftable, cmptable)

char	obsmode[ARB]	# i: Observation mode expression
char	output[ARB]	# i: Output table name
char	photlist[ARB]	# i: List of photometric parameters
real	refwave		# i: Reference wavelength for monchromatic flux
char	wavtable[ARB]	# i: Table containing wavelength array
char	grftable[ARB]	# i: Instrument graph table
char	cmptable[ARB]	# i: Component name table
#--
int	crosslist[MAXPHOT], photvar[MAXPHOT]
int	iobs, nphot, obslen, nwave, degree, iphot, jphot
real	area, phot[MAXPHOT]
pointer	tp, sp, mode, colname, tbname, pcode
pointer	list, par, wave, filt, optr, cptr[MAXPHOT]

data	crosslist  / 11 * 0, 10, 10 /

string	nilstr  ""
string	isempty  "Photometric parameter list is empty"
string	notband  "Expression is not a bandpass"
string	namelist "uresp,pivwv,bandw,fwhm,wpeak,tpeak,avgwv,qtlam,equvw,rectw,\
emflx,refwave,tlambda"

extern	getsynvar
int	lenexpr(), nxtlist(), word_find()
pointer	rdlist(), tbtopn(), locpr(), opn_parmatrix()
real	funit(), pivlam(), rmslam(), fwhmlam(), avglam(), qtlam()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (colname, SZ_FNAME, TY_CHAR)
	call salloc (tbname, SZ_FNAME, TY_CHAR)
	call salloc (pcode, SZ_LINE, TY_INT)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar
	call get_hstarea (area)

	# Read the obsmode string, get its length

	list = rdlist (obsmode)
	obslen = lenexpr (list)

	# Get the photometric parameters to print

	call photnames (namelist, photlist, crosslist, MAXPHOT, photvar, nphot)
	if (nphot == 0)
	    call printerr_str (isempty, obsmode)

	# Initialize the structure that stores the output

	par = opn_parmatrix (list, obslen, namelist, photvar)

	# Get wavelength set 

	call multiwave (wavtable, grftable, cmptable, list, wave, nwave)
	call salloc (filt, nwave, TY_REAL)

	# Create the output table

	if (output[1] != EOS) {
	    tp = tbtopn (output, NEW_FILE, NULL)
	    call tbcdef (tp, optr, "OBSMODE", nilstr, nilstr, -obslen, 1, 1)

	    do iphot = 1, MAXPHOT {
		if (word_find (iphot, namelist, Memc[colname], SZ_FNAME) < 0)
		    break

		if (photvar[iphot] == NO) {
		    cptr[iphot] = NULL

		} else {
		    call strupr (Memc[colname])
		    call tbcdef (tp, cptr[iphot], Memc[colname], 
				 nilstr, nilstr, TY_REAL, 1, 1)
		}
	    }

	    call tbtcre (tp)

	    # Add header parameters

	    call lastfile (grftable, Memc[tbname], SZ_FNAME)
	    call tbhadt (tp, "grftable", Memc[tbname])

	    call lastfile (cmptable, Memc[tbname], SZ_FNAME)
	    call tbhadt (tp, "cmptable", Memc[tbname])

	    call tbhadr (tp, "aperarea", area)
	}

	# Loop over each observation mode expression

	iobs = 1
	while (nxtlist (list, Memc[mode], SZ_LINE) != EOF) {

	    # Compile observation mode into pseudocode

	    call expcompile (Memc[mode], Memi[pcode], SZ_LINE)

	    # Calculate throughput of observation mode

	    call syncalc (Memi[pcode], SZ_COMMAND, locpr(getsynvar), 
			  nwave, Memr[wave], grftable, cmptable, 
			  Memr[filt], degree)

	    if (degree != 0)
		call printerr_str (notband, Memc[mode])

	    # Calculate source independent photometric parameters

	    phot[1] = funit (area, nwave, Memr[wave], Memr[filt])
	    phot[2] = pivlam (nwave, Memr[wave], Memr[filt])
	    phot[3] = rmslam (nwave, Memr[wave], Memr[filt])
	    phot[4] = fwhmlam (nwave, Memr[wave], Memr[filt])

	    call peaklam2 (nwave, Memr[wave], Memr[filt], phot[5], phot[6])

	    phot[7] = avglam (nwave, Memr[wave], Memr[filt])
	    phot[8] = qtlam (nwave, Memr[wave], Memr[filt])

	    call widthlam (nwave, Memr[wave], Memr[filt], 
			   phot[9], phot[10])

	    phot[12] = refwave
	    call monolam (nwave, Memr[wave], Memr[filt], 
			  phot[12], phot[13], phot[11])

	    # Write results to table

	    if (output[1] != EOS)
		call tbeptt (tp, optr, iobs, Memc[mode])

	    jphot = 0
	    do iphot = 1, MAXPHOT {
		if (photvar[iphot] == YES) {
		    jphot = jphot + 1
		    call upd_parmatrix (par, iobs, jphot, phot[iphot])

		    if (output[1] != EOS)
			call tbeptr (tp, cptr[iphot], iobs, phot[iphot])
		}
	    }

	    iobs = iobs + 1
	}

	# Write the results to STDOUT

	call wrt_parmatrix (par)

	# Close files and release memory

	call clssyntab
	if (output[1] != EOS)
	    call tbtclo (tp)

	call mfree (wave, TY_REAL)
	call freelist (list)
	call sfree (sp)
end
