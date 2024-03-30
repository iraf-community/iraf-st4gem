#* HISTORY *
#* B. Simon	23-Sep-94	original

# REPORTGRID -- Report chi squared values for grid spectra

procedure reportgrid (fd, better, best, ngrid, chi2, scale)

int	fd		# i: spool file contain spectrum expressions
int	better		# i: index of second smallest chi squared value
int	best		# i: index of smallest chi squared value
int	ngrid		# i: length of chi2 and scale arrays
real	chi2[ARB]	# i: squared difference of spectrum and data
real	scale[ARB]	# i: scale factor for spectra
#--
int	nc, igrid, irank
pointer	sp, expr
real	relchi2

string	star        "*"
string	header      "       chisq    chisq/best          spectrum\n"
string	footer      "      (*) used in least squares solution\n\n"
string	format      "%2s %11.4g%11.4g   (%s)*%g\n"
string	nospectrum  "reportgrid: spectrum not found, last one read"

int	getline()

begin
	# Allocate temporary string

	call smark (sp)
	call salloc (expr, SZ_LINE, TY_CHAR)

	call seek (fd, BOF)
	call fprintf (STDERR, header)

	do igrid = 1, ngrid  {
	    # Read spectrum expression from spool file, chop off newline

	    nc  = getline (fd, Memc[expr])
	    if (nc <= 0)
		call printerr_str (nospectrum, Memc[expr])

	    Memc[expr+nc-1] = EOS

	    # Compute relative chi squared

	    if (chi2[best] <= 0.0) {
		relchi2 = INDEFR
	    } else {
		relchi2 = chi2[igrid] / chi2[best]
	    }

	    # Display star if spectrum used in fit

	    if (igrid  == better || igrid == best) {
		irank = 1
	    } else {
		irank = 2
	    }

	    # Print results

	    call fprintf (STDERR, format)
	    call pargstr (star[irank])
	    call pargr (chi2[igrid])
	    call pargr (relchi2)
	    call pargstr (Memc[expr])
	    call pargr (scale[igrid])
	}

	call fprintf (STDERR, footer)
	call sfree (sp)
end
