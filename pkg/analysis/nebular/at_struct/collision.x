# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	<tbset.h>
include <math/curfit.h>
include "../at.h"
include "../neberr.h"

define	TE_COL	1

#---------------------------------------------------------------------6 Feb 97--
.help collision.x Apr96 nebular/fivel
.ih
NAME
.nf
 fit_collisions - Fetches collision cross section data and performs fit
       coll_fit - Fit the collision cross section data
       coll_pop - Populate matrix of collision cross sections
  coll_col_find - Fetch column pointers by column name for collision data.   
at_get_func_typ - Get fitting function enum. type from type string.
   at_get_order - Get fit function order from table column units
.fi
.endhelp
#-------------------------------------------------------------------------------
#  FIT_COLLISIONS -- Fetches collision cross section data from table and 
#			performs fit.

procedure fit_collisions (tp, at)

#  Calling arguments:
pointer	tp			# I: atomic data table pointer
pointer	at			# I: atomic data object

#  Declarations:
pointer	colptr			# column pointers
pointer	ft_wt			# array of weights for fit
int	i			# generic
int	ncols			# number of columns in table
int	nrows			# number of rows in table
pointer	nlflag			# array of null flags
int	order			# order of fit
pointer	sp			# top of stack
pointer	temp			# array of temperatures
pointer	coll			# array of collision strengths

string	col_template "Omega(%d->%d)"	# table column name template

#  Functions used:
int 	at_get_order()		# fetch fit order from column attribute
pointer	coll_fit		# perform fit to collision cross sections
int	tbpsta()		# fetch table attributes

#  Memory management:
define	Col_ptr		Memi[colptr+$1-1]
define	Fit_wt		Memr[ft_wt]
define	Nulflg		Memb[nlflag]
define	T_e		Memr[temp+$1-1]
define	Collis		Memr[coll]

errchk	coll_col_find, coll_fit

begin 
	nrows = tbpsta (tp, TBL_NROWS)
	ncols = tbpsta (tp, TBL_NCOLS)

	# Allocate some arrays. 
	call smark (sp)
	call salloc (colptr, ncols, TY_POINTER)

	call salloc (ft_wt,  nrows, TY_REAL)
	call salloc (nlflag, nrows, TY_BOOL)
	call salloc (temp,   nrows, TY_REAL)
	call salloc (coll,   nrows, TY_REAL)

	# Initialize column names. 
	call coll_col_find (tp, col_template, Col_ptr(1), AT_NLVL(at), 
				AT_LOG_TE(at))

	# Initialize the independent array & weights. 
	call tbcgtr (tp, Col_ptr(TE_COL), T_e(1), Nulflg, 1, nrows)
	call amovkr (1., Fit_wt, nrows)

	# Set the limits on T_e.
	AT_TE_MIN(at) = min (T_e(1), T_e(nrows))
	AT_TE_MAX(at) = max (T_e(1), T_e(nrows))
	if (AT_LOG_TE(at)) {
	    AT_TE_MIN(at) = 10. ** (AT_TE_MIN(at))
	    AT_TE_MAX(at) = 10. ** (AT_TE_MAX(at))

	} else {
	    AT_TE_MIN(at) = 1.e4 * AT_TE_MIN(at)
	    AT_TE_MAX(at) = 1.e4 * AT_TE_MAX(at)
	}

	# Perform the fit for each transition. 
	do i = 1, AT_NTRANS(at) {
	    order = at_get_order (Col_ptr(i+TE_COL))
	    if (order < 0)
		order = 0

	    call tbcgtr (tp, Col_ptr(i+TE_COL), Collis, Nulflg, 1, nrows)
	    CURVE(at,i) = coll_fit (T_e(1), Collis, Fit_wt, nrows, 
				AT_TY_FUNC(at), order)
	}

	call sfree (sp)
end


#-------------------------------------------------------------------------------
#  COLL_FIT -- Fit the collision cross section data.

pointer procedure coll_fit (x, y, ft_wt, npts, ty_function, order)

#  Arguments:
real	x[ARB]			# I: independent array (temperatures)
real	y[ARB]			# I: dependent array (collision strengths)
real	ft_wt[ARB]		# I: weights for fit
int	npts			# I: number of points in arrays
int	ty_function		# I: type of fitting function
int	order			# I: order of fitting function
pointer	cv			# O: curve fit structure

#  Local variables:
int	status			# return status of fitting procedure

errchk	cvinit, cvfit

begin
	# Initialize curve fit structure.
	call cvinit (cv, ty_function, order, x[1], x[npts])
	call cvzero (cv)

	# Perform the fit.
	call cvfit (cv, x, y, ft_wt, npts, WTS_UNIFORM, status)
	if (status != OK) 
	    call error (status, "Unable to calculate fit to collision strengths")

	return (cv)
end


#-------------------------------------------------------------------------------
#  COLL_POP -- 	Evaluate fits of collision strength for specified T_e & populate 
#		matrix.  Checks for out-of-bounds T_e

procedure coll_pop (cv, n_trans, coll, n_lvl, t_e)

#  Calling arguments:
pointer	cv[ARB]			# I: curve structures
int	n_trans			# I: no. transitions between atomic levels
double	coll[n_lvl, n_lvl]	# I/O: collision coefficients
int	n_lvl			# I: no. energy levels modelled
real	t_e			# I: electron temperature

#  Local variables:
int	k			# generic
int	lower, upper		# lower, upper level of transition
real	t_min, t_max		# min, max T_e over which curve is defined

#  Functions used:
real	cveval()		# evaluate fitted function
real	cvstatr()		# fetch fit function parameter, TY_REAL

begin
	do upper = 1, nlvl {
	    do lower = 1, nlvl
		coll[upper,lower] = 0.
	}

	k = 0
	do upper = 2, n_lvl {
	    do lower = 1, upper-1 {
		k = k + 1

		# Ensure that evaluation is within T_e bounds for curve.
		t_min = cvstatr (cv[k], CVXMIN)
		t_max = cvstatr (cv[k], CVXMAX)
		if (k <= n_trans) {
		    if (t_e >= t_min && t_e <= t_max) 
		    	coll[upper,lower] = cveval (cv[k], t_e)
		    else
			call error (1, "T_e out of bounds for this ion")
		} else
		    call error (1, "Attempt to use unrecognized transition")
	    }
	}
end


#-------------------------------------------------------------------------------
#  COLL_COL_FIND -- Fetch column pointers by column name for collision data.   

procedure coll_col_find (tp, template, colptr, nlvls, log_te)

#  Calling arguments:
pointer	tp			# I: table descriptor
char	template[ARB]		# I: table column name template
pointer	colptr[ARB]		# O: table column pointers
int	nlvls			# I: no. atomic energy levels expected
bool	log_te			# O: is temperature given in log?

#  Local variables:
char	errmsg[SZ_LINE]		# error message
int	k, lower, upper		# generic
char	name[SZ_COLNAME]	# 

begin
	# Find temperature (independent variable) column.
	log_te = false
	call tbcfnd (tp, "t4", colptr[TE_COL], 1)
	if (colptr[TE_COL] == NULL) {
	    log_te = true
	    call tbcfnd (tp, "log(te)", colptr[TE_COL], 1)
	    if (colptr[TE_COL] == NULL) {
		call sprintf (errmsg, SZ_LINE, 
			"Unable to find temperature column in reference table")
		call error (BAD_REF_DATA, errmsg)
	    }
	}

	# Find collision strength columns by name, or report an error. 
	k = TE_COL
	do upper = 2, nlvls {
	    do lower = 1, upper-1 {
		k = k + 1
	    	call sprintf (name, SZ_COLNAME, template)
		    call pargi (upper)
		    call pargi (lower)
	    	call tbcfnd (tp, name, colptr[k], 1)
	    	if (colptr[k] == NULL) {
		    call sprintf (errmsg, SZ_LINE, 
			"Unable to find column %s in reference table")
		    call pargstr (name)
		    call error (BAD_REF_DATA, errmsg)
		}
	    }
	}

end


#-------------------------------------------------------------------------------
#  AT_GET_FUNC_TYP - Get type of fitting function from type string.

int procedure at_get_func_typ (func_name, max_char)

#  Calling arguments:
char	func_name[ARB]		# I: name of curvefit function
int	max_char		# I: max characters in name
int	func_type		# O: curvfit function type

#  Functions used:
int	strdic()		# string dictionary

begin
	call strlwr (func_name)
	func_type = strdic (func_name, func_name, max_char, CV_FUNCTIONS)

	return (func_type)
end


#-------------------------------------------------------------------------------
#  AT_GET_ORDER - Get fit function order from table column units.

int procedure at_get_order (colptr)

#  Calling arguments:
pointer	colptr			# I: table column pointer
int	order			# O: curvfit order 

#  Local variables:
int	ip, nchar		# generic
char	colunits[SZ_COLUNITS]	# column units

#  Functions used:
int	ctoi()			# convert string to integer value

begin
	call tbcigt (colptr, TBL_COL_UNITS, colunits, SZ_COLUNITS)
	ip = 1
	nchar = ctoi (colunits, ip, order)

	return (order)
end


