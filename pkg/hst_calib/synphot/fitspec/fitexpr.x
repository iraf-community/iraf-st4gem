#* HISTORY *
#* B. Simon	28-Jul-94	original

# FITEXPR -- Compute a least squares fit for a synphot expression

procedure fitexpr (input, expr, output, ftol, maxiter, nprint, 
		   slow, equal, grftable, cmptable, var, nvar)

char	input[ARB]	# i: file containing observed data
char	expr[ARB]	# i: synphot expression used as model
char	output[ARB]	# i: table containing intermediate results
double	ftol		# i: chi squared termination condition
int	maxiter		# i: maximum number of iterations
int	nprint		# i: number of iterations between diagnostic prints
bool	slow		# i: use slow method (simplex) to find lsq fit?
bool	equal		# i: use equal weighting on the data points?
char	grftable[ARB]	# i: instrument graph table
char	cmptable[ARB]	# i: component name table
double	var[ARB]	# u: free variables in least squares fit
int	nvar		# i: number of free variables
#--
int	ndata

extern	fitfunc

begin
	# Initialize synphot table cache and variables

	call inisyntab
	call undefsynvar

	# Set the function to be fit

	call setfitfunc (input, expr, grftable, cmptable, equal, ndata)
			 
	# Compute the values of the variables which minimize 
	# the least squares error of the fit

	if (slow) {
	    call simplexfit (fitfunc, ftol, maxiter, nprint, ndata, nvar, var)
	} else {
	    call levmarfit (fitfunc, ftol, maxiter, nprint, ndata, nvar, var)
	}

	# Write final results and free memory allocated by the fitting process

	call clsfitfunc (output, nvar, var)
	call clssyntab

end
