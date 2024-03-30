define	MAXVAL		4
define	NSPEC		2000

#* HISTORY *
#* B.Simon	18-Jul-96	extracted from calcwave
#* B.Simon	 8-Sep-97	removed int declaration from procedure

# CALCSTEP -- Calculate quadratic coefficients used for wavelength set

procedure calcstep (coef, a, b, c, nwave)

char	coef[ARB]	# i: list of coefficients
real	a		# o: quadratic term
real	b		# o: linear term
real	c		# o: constant term
int	nwave		# o: number of wavelengths
#--
char	lpar, rpar
int	ic, jc, nc, nval
pointer	sp, list, term, value

data	lpar, rpar  / '(' , ')' /

string	badcoef  "Syntax error in coefficient list"

int	stridx(), ctor(), word_fetch()

begin
	call smark (sp)
	call salloc (list, SZ_FNAME, TY_CHAR)
	call salloc (term, SZ_FNAME, TY_CHAR)
	call salloc (value, MAXVAL, TY_REAL)

	# Check to see if coefficient list is enclosed in parentheses
	# Remove them before proceeding

	call strcpy (coef, Memc[list], SZ_FNAME)

	if (Memc[list] != lpar)
	    call printerr_str (badcoef, coef)

	ic = stridx (rpar, Memc[list])

	if (ic == 0)
	    call printerr_str (badcoef, coef)

	Memc[list] = ' '
	Memc[list+ic-1] = ' '

	# Read the coefficient values from the string. The values are:
	# 0: Short wavelength (angstroms)
	# 1: Long wavelength (angstroms)
	# 2: Dispersion at short wavelength (angstroms/pixel)
	# 3: Dispersion at long wavelength (angstroms/pixel)
	# Coefficients 2 and 3 are optional and will be calculated
	# if missing

	ic = 1
	nval = 0
	while (word_fetch (Memc[list], ic, Memc[term], SZ_FNAME) > 0) {
	    if (nval >= MAXVAL)
		call printerr_str (badcoef, coef)

	    jc = 1
	    nc = ctor (Memc[term], jc, Memr[value+nval])

	    if (Memr[value+nval] <= 0.0)
		call printerr_str (badcoef, coef)
		
	    nval = nval + 1
	}

	if (Memr[value] > Memr[value+1])
	    call printerr_str (badcoef, coef)

	# Set any missing coefficients to their default values

	switch (nval) {
	case 2:
	    Memr[value+2] = (Memr[value+1] - Memr[value]) / (NSPEC - 1)
	    Memr[value+3] = Memr[value+2]
	case 3:
	    Memr[value+3] = Memr[value+2]
	case 4:
	    ;
	default:
	    call printerr_str (badcoef, coef)
	}

	# Calculate the number of points in the wavelength set 

	nwave = int (2.0 *  (Memr[value+1] - Memr[value]) / 
		     (Memr[value+3] + Memr[value+2])) + 1

	# Calculate coefficients of quadratic equation used to represent
	# the wavelength set

	c = Memr[value]
	b = Memr[value+2]
	a = (Memr[value+3] * Memr[value+3] - Memr[value+2] * Memr[value+2]) /
	    (4.0 * (Memr[value+1] - Memr[value]))

	call sfree (sp)
end
