#* HISTORY *
#* B. Simon	23-Sep-94	original

# SELECTGRID -- Select spectrum to be fit from grid of spectra

procedure selectgrid (fd, better, best, scale, nvar, var, spectrum, maxch)

int	fd		# i: spool file containing grid spectra
int	better		# i: index of spectrum with second smallest chisq
int	best		# i: index of spectrum with smallest chisq
real	scale[ARB]	# i: scaling factors
int	nvar		# o: number of fit variables
double	var[2]		# o: initial values of fit variables
char	spectrum[ARB]	# o: spectrum expression including fit variables
int	maxch		# i: maximum length of spectrum string
#--
int	nc, igrid
pointer	sp, expr, expr1, expr2

string	notfound  "All spectra are zero in wavelength range of input"
string	specfmt   "$1 * %g * (%s) + (1. - $1) * %g * (%s)"
int	getline()

begin
	# Allocate temporary strings

	call smark (sp)
	call salloc (expr, SZ_LINE, TY_CHAR)
	call salloc (expr1, SZ_LINE, TY_CHAR)
	call salloc (expr2, SZ_LINE, TY_CHAR)

	# Get expressions for better and best spectra

	igrid = 1
	call seek (fd, BOF)

	repeat {
	    nc = getline (fd, Memc[expr])
	    if (nc == EOF)
		break

	    Memc[expr+nc-1] = EOS

	    if (igrid == best) {
		call strcpy (Memc[expr], Memc[expr1], SZ_LINE)
	    } else if (igrid == better) {
		call strcpy (Memc[expr], Memc[expr2], SZ_LINE)
	    }

	    igrid = igrid + 1
	}

	# Create spectrum as sum of better and best spectra

	if (IS_INDEFR (scale[best])) {
	    call error (1, notfound)

	} else if (IS_INDEFR(scale[better])) {
	    nvar = 1
	    var[1] = scale[best]

	    call sprintf (spectrum, maxch, "$1 * (%s)")
	    call pargstr (Memc[expr1])

	} else {
	    nvar = 1
	    var[1] = 0.5 

	    call sprintf (spectrum, maxch, specfmt)
	    call pargr (scale[best])
	    call pargstr (Memc[expr1])
	    call pargr (scale[better])
	    call pargstr (Memc[expr2])
	}

	call sfree (sp)
end

