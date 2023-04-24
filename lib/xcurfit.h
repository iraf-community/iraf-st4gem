# include <c_iraf.h>

void c_dcv_b1leg(double *x, int *order, double *k1, double *k2, double *basis);
/* double	x		# array of data points */
/* int	order		# order of polynomial, order = 1, constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */

void c_dcv_b1cheb(double *x, int *order, double *k1, double *k2, double *basis);
/* double	x		# number of data points */
/* int	order		# order of polynomial, 1 is a constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# array of basis functions */

void c_dcv_b1spline1(double *x, int *npieces, double *k1, double *k2, double *basis, int *left);
/* double	x		# set of data points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */
/* int	left		# index of the appropriate spline functions */

void c_dcv_b1spline3(double *x, int *npieces, double *k1, double *k2, double *basis, int *left);
/* double	x		# array of data points */
/* int	npieces		# number of polynomial pieces */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# array of basis functions */
/* int	left		# array of indices for first non-zero spline */

void c_rcv_b1leg(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# array of data points */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_rcv_b1cheb(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# number of data points */
/* int	order		# order of polynomial, 1 is a constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */

void c_rcv_b1spline1(float *x, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x		# set of data points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */
/* int	left		# index of the appropriate spline functions */

void c_rcv_b1spline3(float *x, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x		# array of data points */
/* int	npieces		# number of polynomial pieces */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */
/* int	left		# array of indices for first non-zero spline */

void c_dcv_bcheb(double *x, int *npts, int *order, double *k1, double *k2, double *basis);
/* double	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */

void c_dcv_bleg(double *x, int *npts, int *order, double *k1, double *k2, double *basis);
/* double	x[npts]		# number of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, 1 is a constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# array of basis functions */

void c_dcv_bspline1(double *x, int *npts, int *npieces, double *k1, double *k2, double *basis, int *left);
/* double	x[npts]		# set of data points */
/* int	npts		# number of points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */
/* int	left[ARB]	# indices of the appropriate spline functions */

void c_dcv_bspline3(double *x, int *npts, int *npieces, double *k1, double *k2, double *basis, int *left);
/* double	x[npts]		# array of data points */
/* int	npts		# number of data points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# array of basis functions */
/* int	left[ARB]	# array of indices for first non-zero spline */

void c_rcv_bcheb(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_rcv_bleg(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# number of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, 1 is a constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */

void c_rcv_bspline1(float *x, int *npts, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x[npts]		# set of data points */
/* int	npts		# number of points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */
/* int	left[ARB]	# indices of the appropriate spline functions */

void c_rcv_bspline3(float *x, int *npts, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x[npts]		# array of data points */
/* int	npts		# number of data points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */
/* int	left[ARB]	# array of indices for first non-zero spline */

void c_dcv_evcheb(double *coeff, double *x, double *yfit, int *npts, int *order, double *k1, double *k2);
/* double	coeff[ARB]		# 1D array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* double	k1, k2			# normalizing constants */

void c_dcv_evleg(double *coeff, double *x, double *yfit, int *npts, int *order, double *k1, double *k2);
/* double	coeff[ARB]		# 1D array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	yfit[npts]		# the fitted points */
/* int	npts			# number of data points */
/* int	order			# order of the polynomial, 1 = constant */
/* double	k1, k2			# normalizing constants */

void c_dcv_evspline1(double *coeff, double *x, double *yfit, int *npts, int *npieces, double *k1, double *k2);
/* double	coeff[ARB]		# array of coefficients */
/* double	x[npts]			# array of x values */
/* double	yfit[npts]		# array of fitted values */
/* int	npts			# number of data points */
/* int	npieces			# number of fitted points minus 1 */
/* double	k1, k2			# normalizing constants */

void c_dcv_evspline3(double *coeff, double *x, double *yfit, int *npts, int *npieces, double *k1, double *k2);
/* double	coeff[ARB]	# array of coeffcients */
/* double	x[npts]		# array of x values */
/* double	yfit[npts]	# array of fitted values */
/* int	npts		# number of data points */
/* int	npieces		# number of polynomial pieces */
/* double	k1, k2		# normalizing constants */

void c_rcv_evcheb(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_rcv_evleg(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of data points */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_rcv_evspline1(float *coeff, float *x, float *yfit, int *npts, int *npieces, float *k1, float *k2);
/* real	coeff[ARB]		# array of coefficients */
/* real	x[npts]			# array of x values */
/* real	yfit[npts]		# array of fitted values */
/* int	npts			# number of data points */
/* int	npieces			# number of fitted points minus 1 */
/* real	k1, k2			# normalizing constants */

void c_rcv_evspline3(float *coeff, float *x, float *yfit, int *npts, int *npieces, float *k1, float *k2);
/* real	coeff[ARB]	# array of coeffcients */
/* real	x[npts]		# array of x values */
/* real	yfit[npts]	# array of fitted values */
/* int	npts		# number of data points */
/* int	npieces		# number of polynomial pieces */
/* real	k1, k2		# normalizing constants */

void c_dcvaccum(IRAFPointer *cv, double *x, double *y, double *w, int *wtflag);
/* pointer	cv		# curve descriptor */
/* double	x		# x value */
/* double	y		# y value */
/* double	w		# weight of the data point */
/* int	wtflag		# type of weighting desired */

void c_cvaccum(IRAFPointer *cv, float *x, float *y, float *w, int *wtflag);
/* pointer	cv		# curve descriptor */
/* real	x		# x value */
/* real	y		# y value */
/* real	w		# weight of the data point */
/* int	wtflag		# type of weighting desired */

void c_dcvacpts(IRAFPointer *cv, double *x, double *y, double *w, int *npts, int *wtflag);
/* pointer	cv		# curve descriptor */
/* double	x[npts]		# array of abcissa */
/* double	y[npts]		# array of ordinates */
/* double	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */

void c_cvacpts(IRAFPointer *cv, float *x, float *y, float *w, int *npts, int *wtflag);
/* pointer	cv		# curve descriptor */
/* real	x[npts]		# array of abcissa */
/* real	y[npts]		# array of ordinates */
/* real	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */

void c_dcvchofac(double *matrix, int *nbands, int *nrows, double *matfac, int *ier);
/* double   matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* double   matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_dcvchoslv(double *matfac, int *nbands, int *nrows, double *vector, double *coeff);
/* double   matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* double   vector[nrows]			# right side of matrix equation */
/* double   coeff[nrows]			# coefficients */

void c_rcvchofac(float *matrix, int *nbands, int *nrows, float *matfac, int *ier);
/* real   matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* real   matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_rcvchoslv(float *matfac, int *nbands, int *nrows, float *vector, float *coeff);
/* real   matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* real   vector[nrows]			# right side of matrix equation */
/* real   coeff[nrows]			# coefficients */

void c_dcvcoeff(IRAFPointer *cv, double *coeff, int *ncoeff);
/* pointer	cv		# curve descriptor */
/* double	coeff[ARB]	# the coefficients of the fit */
/* int	ncoeff		# the number of coefficients */

void c_cvcoeff(IRAFPointer *cv, float *coeff, int *ncoeff);
/* pointer	cv		# curve descriptor */
/* real	coeff[ARB]	# the coefficients of the fit */
/* int	ncoeff		# the number of coefficients */

void c_dcverrors(IRAFPointer *cv, double *y, double *w, double *yfit, int *npts, double *chisqr, double *errors);
/* pointer	cv		# curve descriptor */
/* double	y[ARB]		# data points */
/* double	yfit[ARB]	# fitted data points */
/* double	w[ARB]		# array of weights */
/* int	npts		# number of points */
/* double	chisqr		# reduced chi-squared of fit */
/* double	errors[ARB]	# errors in coefficients */

void c_cverrors(IRAFPointer *cv, float *y, float *w, float *yfit, int *npts, float *chisqr, float *errors);
/* pointer	cv		# curve descriptor */
/* real	y[ARB]		# data points */
/* real	yfit[ARB]	# fitted data points */
/* real	w[ARB]		# array of weights */
/* int	npts		# number of points */
/* real	chisqr		# reduced chi-squared of fit */
/* real	errors[ARB]	# errors in coefficients */

double c_dcveval(IRAFPointer *cv, double *x);
/* pointer	cv		# curve descriptor */
/* double	x		# x value */

float c_cveval(IRAFPointer *cv, float *x);
/* pointer	cv		# curve descriptor */
/* real	x		# x value */

void c_dcvfit(IRAFPointer *cv, double *x, double *y, double *w, int *npts, int *wtflag, int *ier);
/* pointer	cv		# curve descriptor */
/* double	x[npts]		# array of abcissa */
/* double	y[npts]		# array of ordinates */
/* double	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */
/* int	ier		# error code */

void c_cvfit(IRAFPointer *cv, float *x, float *y, float *w, int *npts, int *wtflag, int *ier);
/* pointer	cv		# curve descriptor */
/* real	x[npts]		# array of abcissa */
/* real	y[npts]		# array of ordinates */
/* real	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */
/* int	ier		# error code */

void c_dcvfree(IRAFPointer *cv);
/* pointer	cv	# the curve descriptor */

void c_cvfree(IRAFPointer *cv);
/* pointer	cv	# the curve descriptor */

void c_dcvinit(IRAFPointer *cv, int *curve_type, int *order, double *xmin, double *xmax);
/* pointer	cv		# curve descriptor */
/* int	curve_type	# type of curve to be fitted */
/* int	order		# order of curve to be fitted, or in the case of the */
/* double	xmin		# minimum value of x */
/* double	xmax		# maximum value of x */

void c_cvinit(IRAFPointer *cv, int *curve_type, int *order, float *xmin, float *xmax);
/* pointer	cv		# curve descriptor */
/* int	curve_type	# type of curve to be fitted */
/* int	order		# order of curve to be fitted, or in the case of the */
/* real	xmin		# minimum value of x */
/* real	xmax		# maximum value of x */

void c_dcvpower(IRAFPointer *cv, double *ps_coeff, int *ncoeff);
/* pointer	cv				# Pointer to curfit structure */
/* double	ps_coeff[ncoeff]		# Power series coefficients (output) */
/* int	ncoeff				# Number of coefficients in fit */

void c_dcvepower(IRAFPointer *cv, double *y, double *w, double *yfit, int *npts, double *chisqr, double *perrors);
/* pointer	cv		# curve descriptor */
/* double	y[ARB]		# data points */
/* double	yfit[ARB]	# fitted data points */
/* double	w[ARB]		# array of weights */
/* int	npts		# number of points */
/* double	chisqr		# reduced chi-squared of fit */
/* double	perrors[ARB]	# errors in coefficients */

void c_dcv_mlegen(double *matrix, int *ncoeff);
/* double	matrix[ncoeff, ncoeff] */
/* int	ncoeff */

void c_dcv_etransform(IRAFPointer *cv, double *covar, double *elm, double *perrors, int *ncoeff);
/* pointer	cv */
/* double	covar[ncoeff,ncoeff] */
/* double	elm[ncoeff,ncoeff] */
/* double	perrors[ncoeff] */
/* int	ncoeff */

void c_dcv_legen(double *matrix, double *cf_coeff, double *ps_coeff, int *ncoeff);
/* double	matrix[ncoeff, ncoeff] */
/* double	cf_coeff[ncoeff] */
/* double	ps_coeff[ncoeff] */
/* int	ncoeff */

double c_dcv_legcoeff(int *k, int *n);
/* int	k */
/* int	n */

void c_dcv_mcheby(double *matrix, int *ncoeff);
/* double	matrix[ncoeff, ncoeff]		# Work array for matrix elements */
/* int	ncoeff				# Number of coefficients */

void c_dcv_cheby(double *matrix, double *cf_coeff, double *ps_coeff, int *ncoeff);
/* double	matrix[ncoeff, ncoeff]		# Work array for matrix elements */
/* double	cf_coeff[ncoeff]		# Input curfit coefficients */
/* double	ps_coeff[ncoeff]		# Output power series coefficients */
/* int	ncoeff				# Number of coefficients */

double c_dcv_chebcoeff(int *m, int *n);
/* int	m	# Summation notation index */
/* int	n	# Summation notation index */

void c_dcv_normalize(IRAFPointer *cv, double *ps_coeff, int *ncoeff);
/* pointer	cv			# Pointer to curfit structure */
/* int	ncoeff			# Number of coefficients in fit */
/* double	ps_coeff[ncoeff]	# Power series coefficients */

void c_dcv_enormalize(IRAFPointer *cv, double *elm, int *ncoeff);
/* pointer	cv			# Pointer to curfit structure */
/* double	elm[ncoeff,ncoeff]	# Input transformed matrix */
/* int	ncoeff			# Number of coefficients in fit */

double c_dcv_bcoeff(int *n, int *i);
/* int	n */
/* int	i */

double c_dcv_factorial(int *n);
/* int	n */

void c_cv_mmuld(double *a, double *b, double *c, int *ndim);
/* double   a[ndim,ndim]            #I left input matrix */
/* double   b[ndim,ndim]            #I right input matrix */
/* double   c[ndim,ndim]            #O output matrix */
/* int     ndim                    #I dimensionality of system */

void c_cvpower(IRAFPointer *cv, float *ps_coeff, int *ncoeff);
/* pointer	cv				# Pointer to curfit structure */
/* real	ps_coeff[ncoeff]		# Power series coefficients (output) */
/* int	ncoeff				# Number of coefficients in fit */

void c_cvepower(IRAFPointer *cv, float *y, float *w, float *yfit, int *npts, float *chisqr, float *perrors);
/* pointer	cv		# curve descriptor */
/* real	y[ARB]		# data points */
/* real	yfit[ARB]	# fitted data points */
/* real	w[ARB]		# array of weights */
/* int	npts		# number of points */
/* real	chisqr		# reduced chi-squared of fit */
/* real	perrors[ARB]	# errors in coefficients */

void c_rcv_mlegen(double *matrix, int *ncoeff);
/* double	matrix[ncoeff, ncoeff] */
/* int	ncoeff */

void c_rcv_etransform(IRAFPointer *cv, double *covar, double *elm, float *perrors, int *ncoeff);
/* pointer	cv */
/* double	covar[ncoeff,ncoeff] */
/* double	elm[ncoeff,ncoeff] */
/* real	perrors[ncoeff] */
/* int	ncoeff */

void c_rcv_legen(double *matrix, float *cf_coeff, float *ps_coeff, int *ncoeff);
/* double	matrix[ncoeff, ncoeff] */
/* real	cf_coeff[ncoeff] */
/* real	ps_coeff[ncoeff] */
/* int	ncoeff */

double c_rcv_legcoeff(int *k, int *n);
/* int	k */
/* int	n */

void c_rcv_mcheby(double *matrix, int *ncoeff);
/* double	matrix[ncoeff, ncoeff]		# Work array for matrix elements */
/* int	ncoeff				# Number of coefficients */

void c_rcv_cheby(double *matrix, float *cf_coeff, float *ps_coeff, int *ncoeff);
/* double	matrix[ncoeff, ncoeff]		# Work array for matrix elements */
/* real	cf_coeff[ncoeff]		# Input curfit coefficients */
/* real	ps_coeff[ncoeff]		# Output power series coefficients */
/* int	ncoeff				# Number of coefficients */

double c_rcv_chebcoeff(int *m, int *n);
/* int	m	# Summation notation index */
/* int	n	# Summation notation index */

void c_rcv_normalize(IRAFPointer *cv, float *ps_coeff, int *ncoeff);
/* pointer	cv			# Pointer to curfit structure */
/* int	ncoeff			# Number of coefficients in fit */
/* real	ps_coeff[ncoeff]	# Power series coefficients */

void c_rcv_enormalize(IRAFPointer *cv, double *elm, int *ncoeff);
/* pointer	cv			# Pointer to curfit structure */
/* double	elm[ncoeff,ncoeff]	# Input transformed matrix */
/* int	ncoeff			# Number of coefficients in fit */

double c_rcv_bcoeff(int *n, int *i);
/* int	n */
/* int	i */

double c_rcv_factorial(int *n);
/* int	n */

void c_cv_mmulr(float *a, float *b, float *c, int *ndim);
/* real   a[ndim,ndim]            #I left input matrix */
/* real   b[ndim,ndim]            #I right input matrix */
/* real   c[ndim,ndim]            #O output matrix */
/* int     ndim                    #I dimensionality of system */

void c_dcvrefit(IRAFPointer *cv, double *x, double *y, double *w, int *ier);
/* pointer	cv		# curve descriptor */
/* double	x[ARB]		# x array */
/* double	y[ARB]		# y array */
/* double	w[ARB]		# weight array */
/* int	ier		# error code */

void c_cvrefit(IRAFPointer *cv, float *x, float *y, float *w, int *ier);
/* pointer	cv		# curve descriptor */
/* real	x[ARB]		# x array */
/* real	y[ARB]		# y array */
/* real	w[ARB]		# weight array */
/* int	ier		# error code */

void c_dcvrject(IRAFPointer *cv, double *x, double *y, double *w);
/* pointer	cv		# curve fitting image descriptor */
/* double	x		# x value */
/* double	y		# y value */
/* double	w		# weight of the data point */

void c_cvrject(IRAFPointer *cv, float *x, float *y, float *w);
/* pointer	cv		# curve fitting image descriptor */
/* real	x		# x value */
/* real	y		# y value */
/* real	w		# weight of the data point */

void c_dcvrestore(IRAFPointer *cv, double *fit);
/* pointer	cv		# curve descriptor */
/* double	fit[ARB]	# array containing fit parameters */

void c_cvrestore(IRAFPointer *cv, float *fit);
/* pointer	cv		# curve descriptor */
/* real	fit[ARB]	# array containing fit parameters */

void c_dcvsave(IRAFPointer *cv, double *fit);
/* pointer	cv		# curve descriptor */
/* double	fit[ARB]	# PIXEL array containing curve parameters */

void c_cvsave(IRAFPointer *cv, float *fit);
/* pointer	cv		# curve descriptor */
/* real	fit[ARB]	# PIXEL array containing curve parameters */

void c_dcvset(IRAFPointer *cv, int *curve_type, double *xmin, double *xmax, double *coeff, int *ncoeff);
/* pointer	cv		# curve descriptor */
/* int	curve_type	# the functional form of the curve */
/* double	xmin		# the minimum x value */
/* double	xmax		# the maximum x value */
/* double	coeff[ncoeff]	# the coefficient array */
/* int	ncoeff		# the number of coefficients */

void c_cvset(IRAFPointer *cv, int *curve_type, float *xmin, float *xmax, float *coeff, int *ncoeff);
/* pointer	cv		# curve descriptor */
/* int	curve_type	# the functional form of the curve */
/* real	xmin		# the minimum x value */
/* real	xmax		# the maximum x value */
/* real	coeff[ncoeff]	# the coefficient array */
/* int	ncoeff		# the number of coefficients */

void c_dcvsolve(IRAFPointer *cv, int *ier);
/* pointer	cv 		# curve descriptor */
/* int	ier		# ier = OK, everything OK */

void c_cvsolve(IRAFPointer *cv, int *ier);
/* pointer	cv 		# curve descriptor */
/* int	ier		# ier = OK, everything OK */

int c_dcvstati(IRAFPointer *cv, int *param);
/* pointer	cv			# Curfit pointer */
/* int	param			# Parameter */

double c_dcvstatd(IRAFPointer *cv, int *param);
/* pointer	cv			# Curfit pointer */
/* int	param			# Parameter */

int c_cvstati(IRAFPointer *cv, int *param);
/* pointer	cv			# Curfit pointer */
/* int	param			# Parameter */

float c_cvstatr(IRAFPointer *cv, int *param);
/* pointer	cv			# Curfit pointer */
/* int	param			# Parameter */

void c_dcvvector(IRAFPointer *cv, double *x, double *yfit, int *npts);
/* pointer	cv		# curve descriptor */
/* double	x[npts]		# data x values */
/* double	yfit[npts]	# the fitted y values */
/* int	npts		# number of data points */

void c_cvvector(IRAFPointer *cv, float *x, float *yfit, int *npts);
/* pointer	cv		# curve descriptor */
/* real	x[npts]		# data x values */
/* real	yfit[npts]	# the fitted y values */
/* int	npts		# number of data points */

void c_dcvzero(IRAFPointer *cv);
/* pointer	cv	# pointer to curve descriptor */

void c_cvzero(IRAFPointer *cv);
/* pointer	cv	# pointer to curve descriptor */

