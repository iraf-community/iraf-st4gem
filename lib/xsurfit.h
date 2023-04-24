# include <c_iraf.h>

void c_iscoeff(IRAFPointer *sf, float *coeff, int *ncoeff);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* real	coeff[ARB]	# the coefficients of the fit */
/* int	ncoeff		# the number of coefficients */

float c_iseval(IRAFPointer *sf, float *x, float *y);
/* pointer	sf		# pointer to surface descriptor structure */
/* real	x		# x value */
/* real	y		# y value */

void c_isfree(IRAFPointer *sf);
/* pointer	sf		# pointer to the surface descriptor */

void c_isinit(IRAFPointer *sf, int *surf_type, int *xorder, int *yorder, int *xterms, int *ncols, int *nlines);
/* pointer	sf		# pointer to surface descriptor structure */
/* int	surf_type	# type of surface to be fitted */
/* int	xorder		# x order of surface to be fit, or in the case of the */
/* int	yorder		# y order of surface to be fit, or in the case of the */
/* int	xterms		# cross terms for polynomials? */
/* int	ncols		# number of columns in the surface */
/* int	nlines		# number of lines in the surface */

void c_islaccum(IRAFPointer *sf, int *cols, int *lineno, float *z, float *w, int *ncols, int *wtflag);
/* pointer	sf		# pointer to surface descriptor */
/* int	cols[ncols]	# column values */
/* int	lineno		# lineno of data being accumulated */
/* real	z[ncols]	# surface values on lineno at cols */
/* real	w[ncols]	# weight of the data points */
/* int	ncols		# number of data points */
/* int	wtflag		# type of weighting desired */

void c_islfit(IRAFPointer *sf, int *cols, int *lineno, float *z, float *w, int *ncols, int *wtflag, int *ier);
/* pointer	sf		# pointer to the surface descriptor */
/* int	cols[ncols]	# array of columns */
/* int	lineno    	# lineno */
/* real	z[ncols]	# the surface values */
/* real	w[ncols]	# array of weights */
/* int	ncols		# the number of columns */
/* int	wtflag		# type of weighting */
/* int	ier		# error codes */

void c_islrefit(IRAFPointer *sf, int *cols, int *lineno, float *z, float *w);
/* pointer	sf		# pointer to surface descriptor */
/* int	cols[ARB]	# columns to be fit  */
/* int	lineno		# line number */
/* real	z[ARB]		# surface values */
/* real	w[ARB]		# weight values */

void c_islsolve(IRAFPointer *sf, int *lineno, int *ier);
/* pointer	sf 		# pointer to the surface descriptor structure */
/* int	lineno		# line being fitted in x */
/* int	ier		# ier = 0, everything OK */

void c_islzero(IRAFPointer *sf, int *lineno);
/* pointer	sf	# pointer to the surface descriptor */
/* int	lineno	# line number */

void c_isreplace(IRAFPointer *sf, float *fit);
/* pointer	sf		# surface descriptor */
/* real	fit[ARB]	# array containing the surface parameters and */

void c_isresolve(IRAFPointer *sf, int *lines, int *ier);
/* pointer	sf		# pointer to the surface descriptor structure */
/* int	lines[ARB]	# line numbers included in the fit */
/* int	ier		# error code */

void c_issave(IRAFPointer *sf, float *fit);
/* pointer	sf		# pointer to the surface descriptor */
/* real	fit[ARB]	# array for storing fit */

void c_issolve(IRAFPointer *sf, int *lines, int *nlines, int *ier);
/* pointer	sf		# pointer to the curve descriptor structure */
/* int	lines[ARB]	# line numbers included in the fit */
/* int	nlines		# number of lines fit */
/* int	ier		# error code */

void c_isvector(IRAFPointer *sf, float *x, float *y, float *zfit, int *npts);
/* pointer	sf		# pointer to surface descriptor structure */
/* real	x[ARB]		# x value */
/* real	y[ARB]		# y value */
/* real	zfit[ARB]	# fits surface values */
/* int	npts		# number of data points */

void c_iszero(IRAFPointer *sf);
/* pointer	sf	# pointer to the surface descriptor */

void c_sf_b1leg(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# array of data points */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_sf_b1cheb(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# number of data points */
/* int	order		# order of polynomial, 1 is a constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */

void c_sf_b1spline1(float *x, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x		# set of data points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */
/* int	left		# index of the appropriate spline functions */

void c_sf_b1spline3(float *x, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x		# array of data points */
/* int	npieces		# number of polynomial pieces */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */
/* int	left		# array of indices for first non-zero spline */

void c_sf_bcheb(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_sf_bleg(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# number of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, 1 is a constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */

void c_sf_bspline1(float *x, int *npts, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x[npts]		# set of data points */
/* int	npts		# number of points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */
/* int	left[ARB]	# indices of the appropriate spline functions */

void c_sf_bspline3(float *x, int *npts, int *npieces, float *k1, float *k2, float *basis, int *left);
/* real	x[npts]		# array of data points */
/* int	npts		# number of data points */
/* int	npieces		# number of polynomial pieces minus 1 */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */
/* int	left[ARB]	# array of indices for first non-zero spline */

void c_cv_evcheb(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# EV array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_cv_evleg(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# EV array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of data points */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_cv_evspline1(float *coeff, float *x, float *yfit, int *npts, int *npieces, float *k1, float *k2);
/* real	coeff[ARB]		# array of coefficients */
/* real	x[npts]			# array of x values */
/* real	yfit[npts]		# array of fitted values */
/* int	npts			# number of data points */
/* int	npieces			# number of fitted points minus 1 */
/* real	k1, k2			# normalizing constants */

void c_cv_evspline3(float *coeff, float *x, float *yfit, int *npts, int *npieces, float *k1, float *k2);
/* real	coeff[ARB]	# array of coeffcients */
/* real	x[npts]		# array of x values */
/* real	yfit[npts]	# array of fitted values */
/* int	npts		# number of data points */
/* int	npieces		# number of polynomial pieces */
/* real	k1, k2		# normalizing constants */

void c_sf_evcheb(float *coeff, float *x, float *y, float *zfit, int *npts, int *xterms, int *xorder, int *yorder, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	y[npts] */
/* real	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_sf_evleg(float *coeff, float *x, float *y, float *zfit, int *npts, int *xterms, int *xorder, int *yorder, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	y[npts] */
/* real	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_sf_evspline3(float *coeff, float *x, float *y, float *zfit, int *npts, int *nxpieces, int *nypieces, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# array of coefficients */
/* real	x[npts]			# array of x values */
/* real	y[npts]			# array of y values */
/* real	zfit[npts]		# array of fitted values */
/* int	npts			# number of data points */
/* int	nxpieces, nypieces	# number of fitted points minus 1 */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_sf_evspline1(float *coeff, float *x, float *y, float *zfit, int *npts, int *nxpieces, int *nypieces, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# array of coefficients */
/* real	x[npts]			# array of x values */
/* real	y[npts]			# array of y values */
/* real	zfit[npts]		# array of fitted values */
/* int	npts			# number of data points */
/* int	nxpieces, nypieces	# number of fitted points minus 1 */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_sfchofac(float *matrix, int *nbands, int *nrows, float *matfac, int *ier);
/* real matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* real matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_sfchoslv(float *matfac, int *nbands, int *nrows, float *vector, float *coeff);
/* real matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* real vector[nrows]			# right side of matrix equation */
/* real coeff[nrows]			# coefficients */

