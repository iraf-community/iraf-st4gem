#* HISTORY *
#* B.Simon	30-Sep-94	Original

# TWEENGRID -- Find the best fit between the two best spectra in the grid

procedure tweengrid (spec, nwave, wave, data, weight, output, indegree, units,
		     ftol, maxiter, nprint, slow, grftable, cmptable, 
		     var, nvar)

char 	spec[ARB]	# i: spectrum expression
int	nwave		# i: length of wavelength, data, and weight arrays
real	wave[ARB]	# i: wavelength array
real	data[ARB]	# i: data to be fit
real	weight[ARB]	# i: weights of data paoints
char	output[ARB]	# i: output file name
int	indegree	# i: degree of data
char	units[ARB]	# i: units of data
double	ftol		# i: least squares termination condition
int	maxiter		# i: maximum number of iterations
int	nprint		# i: number of iterations between diagnostic prints
bool	slow		# i: use slow method (simplex) to find lsq fit?
char	grftable[ARB]	# i: instrument graph table
char	cmptable[ARB]	# i: component name table
double	var[ARB]	# u: free variables in least squares fit
int	nvar		# i: number of free variables
#--
extern	gridfit

begin
	# Set the function to be fit

	call setgridfit (spec, grftable, cmptable, nwave, wave, data, 
			  weight, indegree)
			 
	# Compute the values of the variables which minimize 
	# the least squares error of the fit

	if (slow) {
	    call simplexfit (gridfit, ftol, maxiter, nprint, nwave, nvar, var)
	} else {
	    call levmarfit (gridfit, ftol, maxiter, nprint, nwave, nvar, var)
	}

	# Write final results and free memory allocated by the fitting process

	call clsgridfit (output, units, nvar, var, nwave, wave)

end
