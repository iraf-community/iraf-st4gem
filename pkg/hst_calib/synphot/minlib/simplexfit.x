#* HISTORY *
#* B. Simon	12-Sep-94	original

# SIMPLEXFIT -- Nonlinear least squares fit by the downhill simplex method

procedure simplexfit (fitfunc, ftol, maxiter, nprint, ndata, nvar, var)

extern	fitfunc		# i: Function that needs to be fitted
double	ftol		# i: chi squared termination condition
int	maxiter		# i: maximum number of iterations
int	nprint		# i: number of iterations between diagnostic prints
int	ndata		# i: number of data points
int	nvar		# i: number of free variables in fit
double	var[ARB]	# u: free variables
#--
pointer	sp, p, y
real	temp

string	badnvar   "No variables set in fit"
string	badndata  "Insufficient number of data points for fit"
string	badftol   "Fit tolerance is negative"

begin
	# Check input variables

	if (nvar < 1)
	    call printerr_int (badnvar, nvar)

	if (ndata < nvar)
            call printerr_int (badndata, ndata)

        if (ftol < 0.0) {
            temp = ftol
            call printerr_real (badftol, temp)
        }

	# Allocate temporary work arrays

	call smark (sp)
	call salloc (p, nvar*(nvar+1), TY_DOUBLE)
	call salloc (y, nvar+1, TY_DOUBLE)

	# Create initial simplex

	call putsimplex (fitfunc, ndata, nvar, var, Memd[p], Memd[y])

	# Calculate simplex which minimizes residuals

	call simplex (fitfunc, ftol, maxiter, nprint, 
		      ndata, nvar, Memd[p], Memd[y])
	
	# Get final solution from simplex

	call getsimplex (nvar, Memd[p], Memd[y], var)
	call sfree (sp)
end

# GETSIMPLEX -- Copy best result from final simplex to output

procedure getsimplex (nvar, p, y, var)

int	nvar		# i: number of free variables in fit
double	p[nvar+1,nvar]	# i: initial simplex
double	y[nvar+1]	# i: function values at simplex points
double	var[ARB]	# o: free variables
#--
int	i, j, ilo

begin
	# Find smallest value in simplex

	ilo = 1
	do i = 2, nvar+1 {
	    if (y[i] < y[ilo])
		ilo = i
	}

	# Copy vertex with smallest value to output

	do j = 1, nvar
	    var[j] = p[ilo,j]

end

# PUTSIMPLEX -- Build initial simplex for least squares fit

procedure putsimplex (fitfunc, ndata, nvar, var, p, y)

extern	fitfunc		# i: Function that needs to be fitted
int	nvar		# i: number of free variables in fit
int	ndata		# i: number of data points
double	var[ARB]	# i: free variables
double	p[nvar+1,nvar]	# o: initial simplex
double	y[nvar+1]	# o: function values at simplex points
#--
double	step
int	i, j, iflag
pointer	sp, pr, resid

data	step	/ 0.01 /

double	enorm()

begin
	# Allocate temporary arrays

	call smark (sp)
	call salloc (pr, nvar, TY_DOUBLE)
	call salloc (resid, ndata, TY_DOUBLE)

	do i = 1, nvar+1 {
	    # Compute positions of simplex vertices

	    call amovd (var, Memd[pr], nvar)

	    if (i > 1) {
		if (Memd[pr+i-2] == 0.0) {
		    Memd[pr+i-2] = step
		} else {
		    Memd[pr+i-2] = (1.0 + step) * Memd[pr+i-2]
		}
	    }

	    # Compute function values at these vertices

	    iflag = 2
	    call fitfunc (ndata, nvar, Memd[pr], Memd[resid], iflag)
	    y[i] = enorm (ndata, Memd[resid])

	    # Copy vertices into output array

	    do j = 1, nvar
		p[i,j] = Memd[pr+j-1]
	}
	
	call sfree (sp)
end
