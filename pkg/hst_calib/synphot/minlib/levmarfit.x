#* HISTORY *
#* B.Simon	05-Aug-94	Original

# LEVMARFIT -- Nonlinear least squares fit by the Levenberg Marquardt method

procedure levmarfit (fitfunc, xtol, maxiter, nprint, ndata, nvar, var)

extern	fitfunc		# i: Function that needs to be fitted
double	xtol		# i: termination condition
int	maxiter		# i: maximum number of iterations
int	nprint		# i: number of iterations between diagnostic prints
int	ndata		# i: number of data points
int	nvar		# i: number of free variables in fit
double	var[ARB]	# u: free variables
#--
double	ftol, gtol, epsfcn, factor
int	mode, irep, info
pointer	sp, work1, work2, work3, work4, work5, work6
pointer	work7, work8, resid, errmsg
real	temp

string	badnvar   "No variables set in fit"
string	badndata  "Insufficient number of data points for fit"
string	badftol   "Fit tolerance is negative"
string	baditer   "Fit terminated without convergence"

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
        call salloc (work1, nvar, TY_DOUBLE)
        call salloc (work2, nvar*ndata, TY_DOUBLE)
        call salloc (work3, nvar, TY_INT)
        call salloc (work4, nvar, TY_DOUBLE)
        call salloc (work5, nvar, TY_DOUBLE)
        call salloc (work6, nvar, TY_DOUBLE)
        call salloc (work7, nvar, TY_DOUBLE)
        call salloc (work8, ndata, TY_DOUBLE)
        call salloc (resid, ndata, TY_DOUBLE)
        call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Set termination criteria and step sizes for fit

        ftol = 0.0
        gtol = 0.0
        epsfcn = 1.0e-4
        mode = 1
        factor = 0.1

	# Call Levenberg-Marquardt fit routine from minpack

        call lmdif (fitfunc, ndata, nvar, var, Memd[resid], ftol, xtol, gtol,
                    maxiter, epsfcn, Memd[work1], mode, factor, nprint, info,
                    irep, Memd[work2], ndata, Memi[work3], Memd[work4],
                    Memd[work5], Memd[work6], Memd[work7], Memd[work8])

        # Print warning message if convergence not achieved

        if (info == 5) {
            call sprintf (Memc[errmsg], SZ_FNAME, "%d")
            call pargi (maxiter)

            call synphotwarn (baditer, Memc[errmsg])
        }
            
	call sfree (sp)
end
