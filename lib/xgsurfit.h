# include <c_iraf.h>

void c_dgs_b1pol(double *x, int *order, double *k1, double *k2, double *basis);
/* double	x		# data point */
/* int	order		# order of polynomial, order = 1, constant */
/* double	k1, k2		# nomalizing constants, dummy in this case */
/* double	basis[ARB]	# basis functions */

void c_dgs_b1leg(double *x, int *order, double *k1, double *k2, double *basis);
/* double	x		# data point */
/* int	order		# order of polynomial, order = 1, constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */

void c_dgs_b1cheb(double *x, int *order, double *k1, double *k2, double *basis);
/* double	x		# number of data points */
/* int	order		# order of polynomial, 1 is a constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# array of basis functions */

void c_rgs_b1pol(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# data point */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# nomalizing constants, dummy in this case */
/* real	basis[ARB]	# basis functions */

void c_rgs_b1leg(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# data point */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_rgs_b1cheb(float *x, int *order, float *k1, float *k2, float *basis);
/* real	x		# number of data points */
/* int	order		# order of polynomial, 1 is a constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */

void c_dgs_bpol(double *x, int *npts, int *order, double *k1, double *k2, double *basis);
/* double	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */

void c_dgs_bcheb(double *x, int *npts, int *order, double *k1, double *k2, double *basis);
/* double	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# basis functions */

void c_dgs_bleg(double *x, int *npts, int *order, double *k1, double *k2, double *basis);
/* double	x[npts]		# number of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, 1 is a constant */
/* double	k1, k2		# normalizing constants */
/* double	basis[ARB]	# array of basis functions */

void c_rgs_bpol(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_rgs_bcheb(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# array of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, order = 1, constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# basis functions */

void c_rgs_bleg(float *x, int *npts, int *order, float *k1, float *k2, float *basis);
/* real	x[npts]		# number of data points */
/* int	npts		# number of points */
/* int	order		# order of polynomial, 1 is a constant */
/* real	k1, k2		# normalizing constants */
/* real	basis[ARB]	# array of basis functions */

void c_dgschofac(double *matrix, int *nbands, int *nrows, double *matfac, int *ier);
/* double matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* double matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_dgschoslv(double *matfac, int *nbands, int *nrows, double *vector, double *coeff);
/* double matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* double vector[nrows]			# right side of matrix equation */
/* double coeff[nrows]			# coefficients */

void c_rgschofac(float *matrix, int *nbands, int *nrows, float *matfac, int *ier);
/* real matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* real matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_rgschoslv(float *matfac, int *nbands, int *nrows, float *vector, float *coeff);
/* real matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* real vector[nrows]			# right side of matrix equation */
/* real coeff[nrows]			# coefficients */

void c_dgs_1devpoly(double *coeff, double *x, double *yfit, int *npts, int *order, double *k1, double *k2);
/* double	coeff[ARB]		# EV array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* double	k1, k2			# normalizing constants */

void c_dgs_1devcheb(double *coeff, double *x, double *yfit, int *npts, int *order, double *k1, double *k2);
/* double	coeff[ARB]		# EV array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* double	k1, k2			# normalizing constants */

void c_dgs_1devleg(double *coeff, double *x, double *yfit, int *npts, int *order, double *k1, double *k2);
/* double	coeff[ARB]		# EV array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	yfit[npts]		# the fitted points */
/* int	npts			# number of data points */
/* int	order			# order of the polynomial, 1 = constant */
/* double	k1, k2			# normalizing constants */

void c_rgs_1devpoly(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# EV array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_rgs_1devcheb(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# EV array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_rgs_1devleg(float *coeff, float *x, float *yfit, int *npts, int *order, float *k1, float *k2);
/* real	coeff[ARB]		# EV array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	yfit[npts]		# the fitted points */
/* int	npts			# number of data points */
/* int	order			# order of the polynomial, 1 = constant */
/* real	k1, k2			# normalizing constants */

void c_dgs_evpoly(double *coeff, double *x, double *y, double *zfit, int *npts, int *xterms, int *xorder, int *yorder, double *k1x, double *k2x, double *k1y, double *k2y);
/* double	coeff[ARB]		# 1D array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	y[npts] */
/* double	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* double	k1x, k2x		# normalizing constants */
/* double	k1y, k2y */

void c_dgs_evcheb(double *coeff, double *x, double *y, double *zfit, int *npts, int *xterms, int *xorder, int *yorder, double *k1x, double *k2x, double *k1y, double *k2y);
/* double	coeff[ARB]		# 1D array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	y[npts] */
/* double	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* double	k1x, k2x		# normalizing constants */
/* double	k1y, k2y */

void c_dgs_evleg(double *coeff, double *x, double *y, double *zfit, int *npts, int *xterms, int *xorder, int *yorder, double *k1x, double *k2x, double *k1y, double *k2y);
/* double	coeff[ARB]		# 1D array of coefficients */
/* double	x[npts]			# x values of points to be evaluated */
/* double	y[npts] */
/* double	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* double	k1x, k2x		# normalizing constants */
/* double	k1y, k2y */

void c_gs_asumvpd(double *a, double *b, double *c, double *d, int *npts);
/* double	a[ARB]		# first input vector */
/* double	b[ARB]		# second input vector */
/* double	c[ARB]		# third vector */
/* double	d[ARB]		# output vector */
/* int	npts		# number of points */

void c_rgs_evpoly(float *coeff, float *x, float *y, float *zfit, int *npts, int *xterms, int *xorder, int *yorder, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	y[npts] */
/* real	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_rgs_evcheb(float *coeff, float *x, float *y, float *zfit, int *npts, int *xterms, int *xorder, int *yorder, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	y[npts] */
/* real	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_rgs_evleg(float *coeff, float *x, float *y, float *zfit, int *npts, int *xterms, int *xorder, int *yorder, float *k1x, float *k2x, float *k1y, float *k2y);
/* real	coeff[ARB]		# 1D array of coefficients */
/* real	x[npts]			# x values of points to be evaluated */
/* real	y[npts] */
/* real	zfit[npts]		# the fitted points */
/* int	npts			# number of points to be evaluated */
/* int	xterms			# cross terms ? */
/* int	xorder,yorder		# order of the polynomials in x and y */
/* real	k1x, k2x		# normalizing constants */
/* real	k1y, k2y */

void c_gs_asumvpr(float *a, float *b, float *c, float *d, int *npts);
/* real	a[ARB]		# first input vector */
/* real	b[ARB]		# second input vector */
/* real	c[ARB]		# third vector */
/* real	d[ARB]		# output vector */
/* int	npts		# number of points */

void c_dgsaccum(IRAFPointer *sf, double *x, double *y, double *z, double *w, int *wtflag);
/* pointer	sf		# surface descriptor */
/* double	x		# x value */
/* double	y		# y value */
/* double	z		# z value */
/* double	w		# weight */
/* int	wtflag		# type of weighting */

void c_gsaccum(IRAFPointer *sf, float *x, float *y, float *z, float *w, int *wtflag);
/* pointer	sf		# surface descriptor */
/* real	x		# x value */
/* real	y		# y value */
/* real	z		# z value */
/* real	w		# weight */
/* int	wtflag		# type of weighting */

void c_dgsacpts(IRAFPointer *sf, double *x, double *y, double *z, double *w, int *npts, int *wtflag);
/* pointer	sf		# surface descriptor */
/* double	x[npts]		# array of x values */
/* double	y[npts]		# array of y values */
/* double	z[npts]		# data array */
/* double	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */

void c_gsacpts(IRAFPointer *sf, float *x, float *y, float *z, float *w, int *npts, int *wtflag);
/* pointer	sf		# surface descriptor */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	z[npts]		# data array */
/* real	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */

void c_dgsadd(IRAFPointer *sf1, IRAFPointer *sf2, IRAFPointer *sf3);
/* pointer	sf1		# pointer to the first surface */
/* pointer	sf2		# pointer to the second surface */
/* pointer	sf3		# pointer to the output surface */

void c_gsadd(IRAFPointer *sf1, IRAFPointer *sf2, IRAFPointer *sf3);
/* pointer	sf1		# pointer to the first surface */
/* pointer	sf2		# pointer to the second surface */
/* pointer	sf3		# pointer to the output surface */

void c_dgscoeff(IRAFPointer *sf, double *coeff, int *ncoeff);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* double	coeff[ARB]	# the coefficients of the fit */
/* int	ncoeff		# the number of coefficients */

void c_gscoeff(IRAFPointer *sf, float *coeff, int *ncoeff);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* real	coeff[ARB]	# the coefficients of the fit */
/* int	ncoeff		# the number of coefficients */

void c_dgscopy(IRAFPointer *sf1, IRAFPointer *sf2);
/* pointer	sf1		# pointer to original surface */
/* pointer	sf2		# pointer to the new surface */

void c_gscopy(IRAFPointer *sf1, IRAFPointer *sf2);
/* pointer	sf1		# pointer to original surface */
/* pointer	sf2		# pointer to the new surface */

void c_dgsder(IRAFPointer *sf1, double *x, double *y, double *zfit, int *npts, int *nxd, int *nyd);
/* pointer	sf1		# pointer to the previous surface */
/* double	x[npts]		# x values */
/* double	y[npts]		# y values */
/* double	zfit[npts]	# fitted values */
/* int	npts		# number of points */
/* int	nxd, nyd	# order of the derivatives in x and y */

void c_gsder(IRAFPointer *sf1, float *x, float *y, float *zfit, int *npts, int *nxd, int *nyd);
/* pointer	sf1		# pointer to the previous surface */
/* real	x[npts]		# x values */
/* real	y[npts]		# y values */
/* real	zfit[npts]	# fitted values */
/* int	npts		# number of points */
/* int	nxd, nyd	# order of the derivatives in x and y */

void c_dgserrors(IRAFPointer *sf, double *z, double *w, double *zfit, double *chisqr, double *errors);
/* pointer	sf		# curve descriptor */
/* double	z[ARB]		# data points */
/* double	w[ARB]		# array of weights */
/* double	zfit[ARB]	# fitted data points */
/* double	chisqr		# reduced chi-squared of fit */
/* double	errors[ARB]	# errors in coefficients */

void c_gserrors(IRAFPointer *sf, float *z, float *w, float *zfit, float *chisqr, float *errors);
/* pointer	sf		# curve descriptor */
/* real	z[ARB]		# data points */
/* real	w[ARB]		# array of weights */
/* real	zfit[ARB]	# fitted data points */
/* real	chisqr		# reduced chi-squared of fit */
/* real	errors[ARB]	# errors in coefficients */

double c_dgseval(IRAFPointer *sf, double *x, double *y);
/* pointer	sf		# pointer to surface descriptor structure */
/* double	x		# x value */
/* double	y		# y value */

float c_gseval(IRAFPointer *sf, float *x, float *y);
/* pointer	sf		# pointer to surface descriptor structure */
/* real	x		# x value */
/* real	y		# y value */

void c_dgsfit(IRAFPointer *sf, double *x, double *y, double *z, double *w, int *npts, int *wtflag, int *ier);
/* pointer	sf		# surface descriptor */
/* double	x[npts]		# array of x values */
/* double	y[npts]		# array of y values */
/* double	z[npts]		# data array */
/* double	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */
/* int	ier		# ier = OK, everything OK */

void c_gsfit(IRAFPointer *sf, float *x, float *y, float *z, float *w, int *npts, int *wtflag, int *ier);
/* pointer	sf		# surface descriptor */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	z[npts]		# data array */
/* real	w[npts]		# array of weights */
/* int	npts		# number of data points */
/* int	wtflag		# type of weighting */
/* int	ier		# ier = OK, everything OK */

void c_dgsfree(IRAFPointer *sf);
/* pointer	sf	# the surface descriptor */

void c_gsfree(IRAFPointer *sf);
/* pointer	sf	# the surface descriptor */

double c_dgsgcoeff(IRAFPointer *sf, int *xorder, int *yorder);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* int	xorder		# X order of desired coefficent */
/* int	yorder		# Y order of desired coefficent */

float c_gsgcoeff(IRAFPointer *sf, int *xorder, int *yorder);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* int	xorder		# X order of desired coefficent */
/* int	yorder		# Y order of desired coefficent */

void c_dgsinit(IRAFPointer *sf, int *surface_type, int *xorder, int *yorder, int *xterms, double *xmin, double *xmax, double *ymin, double *ymax);
/* pointer	sf		# surface descriptor */
/* int	surface_type	# type of surface to be fitted */
/* int	xorder		# x order of surface to be fit */
/* int	yorder		# y order of surface to be fit */
/* int	xterms		# presence of cross terms */
/* double   xmin		# minimum value of x */
/* double	xmax		# maximum value of x */
/* double	ymin		# minimum value of y */
/* double	ymax		# maximum value of y */

void c_gsinit(IRAFPointer *sf, int *surface_type, int *xorder, int *yorder, int *xterms, float *xmin, float *xmax, float *ymin, float *ymax);
/* pointer	sf		# surface descriptor */
/* int	surface_type	# type of surface to be fitted */
/* int	xorder		# x order of surface to be fit */
/* int	yorder		# y order of surface to be fit */
/* int	xterms		# presence of cross terms */
/* real   xmin		# minimum value of x */
/* real	xmax		# maximum value of x */
/* real	ymin		# minimum value of y */
/* real	ymax		# maximum value of y */

void c_dgsrefit(IRAFPointer *sf, double *x, double *y, double *z, double *w, int *ier);
/* pointer	sf		# surface descriptor */
/* double	x[ARB]		# array of x values */
/* double	y[ARB]		# array of y values */
/* double	z[ARB]		# data array */
/* double	w[ARB]		# array of weights */
/* int	ier		# ier = OK, everything OK */

void c_gsrefit(IRAFPointer *sf, float *x, float *y, float *z, float *w, int *ier);
/* pointer	sf		# surface descriptor */
/* real	x[ARB]		# array of x values */
/* real	y[ARB]		# array of y values */
/* real	z[ARB]		# data array */
/* real	w[ARB]		# array of weights */
/* int	ier		# ier = OK, everything OK */

void c_dgsrej(IRAFPointer *sf, double *x, double *y, double *z, double *w, int *wtflag);
/* pointer	sf		# surface descriptor */
/* double	x		# x value */
/* double	y		# y value */
/* double	z		# z value */
/* double	w		# weight */
/* int	wtflag		# type of weighting */

void c_gsrej(IRAFPointer *sf, float *x, float *y, float *z, float *w, int *wtflag);
/* pointer	sf		# surface descriptor */
/* real	x		# x value */
/* real	y		# y value */
/* real	z		# z value */
/* real	w		# weight */
/* int	wtflag		# type of weighting */

void c_dgsrestore(IRAFPointer *sf, double *fit);
/* pointer	sf		# surface descriptor */
/* double	fit[ARB]	# array containing the surface parameters and */

void c_gsrestore(IRAFPointer *sf, float *fit);
/* pointer	sf		# surface descriptor */
/* real	fit[ARB]	# array containing the surface parameters and */

void c_dgssave(IRAFPointer *sf, double *fit);
/* pointer	sf		# pointer to the surface descriptor */
/* double	fit[ARB]	# array for storing fit */

void c_gssave(IRAFPointer *sf, float *fit);
/* pointer	sf		# pointer to the surface descriptor */
/* real	fit[ARB]	# array for storing fit */

void c_dgsscoeff(IRAFPointer *sf, int *xorder, int *yorder, double *coeff);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* int	xorder		# X order of desired coefficent */
/* int	yorder		# Y order of desired coefficent */
/* double	coeff		# Coefficient value */

void c_gsscoeff(IRAFPointer *sf, int *xorder, int *yorder, float *coeff);
/* pointer	sf		# pointer to the surface fitting descriptor */
/* int	xorder		# X order of desired coefficent */
/* int	yorder		# Y order of desired coefficent */
/* real	coeff		# Coefficient value */

void c_dgssolve(IRAFPointer *sf, int *ier);
/* pointer	sf 		# curve descriptor */
/* int	ier		# ier = OK, everything OK */

void c_gssolve(IRAFPointer *sf, int *ier);
/* pointer	sf 		# curve descriptor */
/* int	ier		# ier = OK, everything OK */

double c_dgsgetd(IRAFPointer *sf, int *parameter);
/* pointer	sf		# pointer to the surface fit */
/* int	parameter	# parameter to be fetched */

int c_dgsgeti(IRAFPointer *sf, int *parameter);
/* pointer sf		# pointer to the surface fit */
/* int	parameter	# integer parameter */

float c_gsgetr(IRAFPointer *sf, int *parameter);
/* pointer	sf		# pointer to the surface fit */
/* int	parameter	# parameter to be fetched */

int c_gsgeti(IRAFPointer *sf, int *parameter);
/* pointer sf		# pointer to the surface fit */
/* int	parameter	# integer parameter */

void c_dgssub(IRAFPointer *sf1, IRAFPointer *sf2, IRAFPointer *sf3);
/* pointer	sf1		# pointer to the first surface */
/* pointer	sf2		# pointer to the second surface */
/* pointer	sf3		# pointer to the output surface */

void c_gssub(IRAFPointer *sf1, IRAFPointer *sf2, IRAFPointer *sf3);
/* pointer	sf1		# pointer to the first surface */
/* pointer	sf2		# pointer to the second surface */
/* pointer	sf3		# pointer to the output surface */

void c_dgsvector(IRAFPointer *sf, double *x, double *y, double *zfit, int *npts);
/* pointer	sf		# pointer to surface descriptor structure */
/* double	x[ARB]		# x value */
/* double	y[ARB]		# y value */
/* double	zfit[ARB]	# fits surface values */
/* int	npts		# number of data points */

void c_gsvector(IRAFPointer *sf, float *x, float *y, float *zfit, int *npts);
/* pointer	sf		# pointer to surface descriptor structure */
/* real	x[ARB]		# x value */
/* real	y[ARB]		# y value */
/* real	zfit[ARB]	# fits surface values */
/* int	npts		# number of data points */

void c_dgszero(IRAFPointer *sf);
/* pointer	sf	# pointer to surface descriptor */

void c_gszero(IRAFPointer *sf);
/* pointer	sf	# pointer to surface descriptor */

