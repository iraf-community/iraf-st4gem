# include <c_iraf.h>

void c_arbpix(float *datain, float *dataout, int *npts, int *interp_type, int *boundary_type);
/* real	datain[ARB]		# input data array */
/* real	dataout[ARB]		# output data array - cannot be same as datain */
/* int	npts			# no. of data points */
/* int	interp_type		# interpolator type */
/* int	boundary_type		# boundary type, at present must be BOUNDARY_EXT */

float c_ii_badpix(float *datain, int *npix, int *index, int *interp_type);
/* real	datain[ARB]	# datain array, y[1] and y[n] guaranteed to be good */
/* int	npix		# length of y array */
/* int	index		# index of bad point to replace */
/* int	interp_type	# interpolator type */

float c_ii_newpix(float *x, float *xarray, float *data, int *npts, int *index, int *interp_type);
/* real	x		# point to interpolate */
/* real	xarray[ARB]	# x values */
/* real	data[ARB]	# data values */
/* int 	npts		# size of data array */
/* int	index		# index such that xarray[index] < x < xarray[index+1] */
/* int 	interp_type	# interpolator type */

void c_ii_badsinc(float *datain, float *dataout, int *npts, int *nsinc, int *ntaper, float *staper, float *min_bdx);
/* real	datain[ARB]	# input data including bad pixels with INDEF values */
/* real	dataout[ARB]	# output data  */
/* int	npts		# number of data values */
/* int	nsinc		# sinc truncation length */
/* int	ntaper		# start of triangular taper */
/* real	staper		# slope of triangular taper */
/* real	min_bdx		# minimum  distance from interpolation point */

void c_arider(float *x, float *datain, int *npix, float *derivs, int *nder, int *interp_type);
/* real	x		# need 1 <= x <= n */
/* real	datain[ARB]	# data values */
/* int	npix		# number of data values */
/* real	derivs[ARB]	# derivatives out -- derivs[1] is function value */
/* int	nder		# total number of values returned in derivs */
/* int	interp_type	# type of interpolator */

float c_arieval(float *x, float *datain, int *npts, int *interp_type);
/* real	x		# x value, 1 <= x <= n */
/* real	datain[ARB]	# array of data values */
/* int	npts		# number of data values */
/* int	interp_type	# interpolant type */

void c_asider(IRAFPointer *asi, float *x, float *der, int *nder);
/* pointer	asi		# interpolant descriptor */
/* real	x		# x value */
/* real	der[ARB]	# derivatives, der[1] is value der[2] is f prime  */
/* int	nder		# number items returned = 1 + number of derivatives */

float c_asieval(IRAFPointer *asi, float *x);
/* pointer	asi		# interpolator descriptor */
/* real x			# x value */

void c_asifit(IRAFPointer *asi, float *datain, int *npix);
/* pointer	asi		# interpolant descriptor */
/* real	datain[ARB]	# data array */
/* int	npix		# nunber of data points */

void c_asifree(IRAFPointer *asi);
/* pointer	asi	# interpolant descriptor */

int c_asigeti(IRAFPointer *asi, int *param);
/* pointer	asi		# interpolant descriptor */
/* int	param		# parameter to be fetched */

float c_asigrl(IRAFPointer *asi, float *a, float *b);
/* pointer	asi		# interpolant descriptor */
/* real	a		# lower limit for integral */
/* real	b		# upper limit for integral */

void c_asiinit(IRAFPointer *asi, int *interp_type);
/* pointer	asi		# interpolant descriptor */
/* int	interp_type	# interpolant type */

void c_asirestore(IRAFPointer *asi, float *interpolant);
/* pointer	asi			# interpolant descriptor */
/* real	interpolant[ARB]	# array containing the interpolant */

void c_asisave(IRAFPointer *asi, float *interpolant);
/* pointer	asi			# interpolant descriptor */
/* real	interpolant[ARB]	# array containing the interpolant */

void c_asivector(IRAFPointer *asi, float *x, float *y, int *npix);
/* pointer	asi		# interpolator descriptor */
/* real	x[ARB]		# ordered x array */
/* real	y[ARB]		# interpolated values */
/* int	npix		# no. of points in x */

float c_ii_1dinteg(float *coeff, float *a, float *b, int *interp_type);
/* real	coeff[ARB]	# 1D array of coefficients */
/* real	a		# lower limit for integral */
/* real	b		# upper limit for integral */
/* int	interp_type	# type of 1D interpolant */

void c_ii_getpcoeff(float *coeff, int *index, float *pcoeff, int *interp_type);
/* real	coeff[ARB]	# coefficient array */
/* int	index		# coefficients wanted for index < x < index + 1 */
/* real	pcoeff[ARB]	# polynomial coefficients */
/* int	interp_type	# type of interpolant */

void c_ii_sincigrl(float *a, float *b, float *sum, float *data, int *npix, int *nsinc, int *ntaper, float *staper, float *mindx);
/* real	a, b			# integral limits */
/* real	sum			# output integral value */
/* real	data[npix]		# input data array */
/* int	npix			# number of pixels */
/* int	nsinc			# sinc truncation length */
/* int	ntaper			# start of triangular taper */
/* real	staper			# slope of triangular taper */
/* real	mindx			# interpolation minimum */

void c_ii_binearest(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *npts);
/* real	coeff[ARB]	# 1D coefficient array */
/* int	first_point	# offset of first data point */
/* int	len_coeff	# row length of coeff */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	zfit[npts]	# array of interpolated values */
/* int	npts		# number of points to be evaluated */

void c_ii_bilinear(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *npts);
/* real	coeff[ARB]	# 1D array of coefficients */
/* int	first_point	# offset of first data point */
/* int	len_coeff	# row length of coeff */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	zfit[npts]	# array of interpolated values */
/* int	npts		# number of data points */

void c_ii_bipoly3(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *npts);
/* real	coeff[ARB]	# 1D array of coefficients */
/* int	first_point	# offset first point */
/* int	len_coeff	# row length of the coefficient array */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	zfit[npts]	# array of fitted values */
/* int	npts		# number of points to be evaluated */

void c_ii_bipoly5(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *npts);
/* real	coeff[ARB]	# 1D array of coefficients */
/* int	first_point	# offset to first data point */
/* int	len_coeff	# row length of coeff */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	zfit[npts]	# array of fitted values */
/* int	npts		# number of points */

void c_ii_bispline3(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *npts);
/* real	coeff[ARB]	# 1D array of coefficients */
/* int	first_point	# offset to first data point */
/* int	len_coeff	# row length of coeff */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	zfit[npts]	# array of interpolated values */
/* int	npts		# number of points to be evaluated */

void c_ii_nearest(float *x, float *y, int *npts, float *data);
/* real	x[ARB]		# x values, must be within [1,npts] */
/* real	y[ARB]		# interpolated values returned to user */
/* int	npts		# number of x values */
/* real	data[ARB]	# data to be interpolated */

void c_ii_linear(float *x, float *y, int *npts, float *data);
/* real	x[ARB]		# x values, must be within [1,npts] */
/* real	y[ARB]		# interpolated values returned to user */
/* int	npts		# number of x values */
/* real	data[ARB]	# data to be interpolated */

void c_ii_poly3(float *x, float *y, int *npts, float *data);
/* real	x[ARB]		# x values, must be within [1,npts] */
/* real	y[ARB]		# interpolated values returned to user */
/* int	npts		# number of x values */
/* real	data[ARB]	# data to be interpolated from a[0] to a[npts+2] */

void c_ii_poly5(float *x, float *y, int *npts, float *data);
/* real	x[ARB]		# x values, must be within [1,npts] */
/* real	y[ARB]		# interpolated values returned to user */
/* int	npts		# number of x values */
/* real	data[ARB]	# data to be interpolated - from a[-1] to a[npts+3] */

void c_ii_spline3(float *x, float *y, int *npts, float *bcoeff);
/* real	x[ARB]		# x values, must be within [1,npts] */
/* real	y[ARB]		# interpolated values returned to user */
/* int	npts		# number of x values */
/* real	bcoeff[ARB]	# basis spline coefficients - from a[0] to a[npts+1] */

void c_ii_sinc(float *x, float *y, int *npts, float *data, int *npix, int *nsinc, int *ntaper, float *staper, float *mindx);
/* real	x[ARB]		# x values, must be within [1,npts] */
/* real	y[ARB]		# interpolated values returned to user */
/* int	npts		# number of x values */
/* real	data[ARB]	# data to be interpolated */
/* int	npix		# number of data pixels */
/* int	nsinc		# sinc truncation length */
/* int	ntaper		# start of triangular taper */
/* real	staper		# slope of triangular taper */
/* real	mindx		# interpolation minimum */

void c_ii_grnearest(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *nxpts, int *nypts, int *len_zfit);
/* real	coeff[ARB]		# 1D coefficient array */
/* int	first_point		# offset of first data point */
/* int	len_coeff		# row length of coeff  */
/* real	x[nxpts]		# array of x values */
/* real	y[nypts]		# array of y values */
/* real	zfit[len_zfit,ARB]	# array of interpolatedvalues */
/* int	nxpts			# number of x values */
/* int	nypts			# number of y values */
/* int	len_zfit		# row length of zfit */

void c_ii_grlinear(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *nxpts, int *nypts, int *len_zfit);
/* real	coeff[ARB]		# 1D array of coefficients */
/* int	first_point		# offset of first data point */
/* int	len_coeff		# row length of coeff */
/* real	x[nxpts]		# array of x values */
/* real	y[nypts]		# array of y values */
/* real	zfit[len_zfit,ARB]	# array of interpolated values */
/* int	nxpts			# number of x values */
/* int	nypts			# number of y values */
/* int	len_zfit		# row length of zfit */

void c_ii_grpoly3(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *nxpts, int *nypts, int *len_zfit);
/* real	coeff[ARB]		# 1D array of coefficients */
/* int	first_point		# offset of first data point */
/* int	len_coeff		# length of row of coeffcient */
/* real	x[nxpts]		# array of x values */
/* real	y[nypts]		# array of y values */
/* real	zfit[len_zfit,ARB]	# array of interpolatedvalues */
/* int	nxpts			# number of x points */
/* int	nypts			# number of y points */
/* int	len_zfit		# row length of zfit */

void c_ii_grpoly5(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *nxpts, int *nypts, int *len_zfit);
/* real	coeff[ARB]		# 1D array of coefficients */
/* int	first_point		# offset of first data point */
/* int	len_coeff		# row length of coeff */
/* real	x[nxpts]		# array of x values */
/* real	y[nypts]		# array of y values */
/* real	zfit[len_zfit,ARB]	# array of fitted values */
/* int	nxpts			# number of x points */
/* int	nypts			# number of y points */
/* int	len_zfit		# row length of zfit */

void c_ii_grspline3(float *coeff, int *first_point, int *len_coeff, float *x, float *y, float *zfit, int *nxpts, int *nypts, int *len_zfit);
/* real	coeff[ARB]		# 1D array of coefficients */
/* int	first_point		# offset of first data point */
/* int	len_coeff		# row length of coeff */
/* real	x[nxpts]		# array of x values */
/* real	y[nypts]		# array of y values */
/* real	zfit[len_zfit,ARB]	# array of interpolated values */
/* int	nxpts			# number of x values */
/* int	nypts			# number of y values */
/* int	len_zfit		# row length of zfit */

void c_ia_pcpoly3(float *x, float *datain, int *npts, float *pcoeff);
/* real	x		# x value */
/* real	datain[ARB]	# array of input data */
/* int	npts		# number of data points */
/* real	pcoeff[ARB]	# array of polynomial coefficients */

void c_ia_pcpoly5(float *x, float *datain, int *npts, float *pcoeff);
/* real	x		# x value */
/* real	datain[ARB]	# array of input data */
/* int	npts		# number of data points */
/* real	pcoeff[ARB]	# array of polynomial coefficients */

void c_ia_pcspline3(float *x, float *datain, int *npts, float *pcoeff);
/* real	x		# x value */
/* real	datain[ARB]	# data array */
/* int	npts		# number of data points */
/* real	pcoeff[ARB]	# array of polynomial coefficients */

void c_ii_sincder(float *x, float *der, int *nder, float *data, int *npix, int *nsinc, int *ntaper, float *staper, float *mindx);
/* real	x		# x value */
/* real	der[ARB]	# derivatives to return */
/* int	nder		# number of derivatives */
/* real	data[npix]	# data to be interpolated */
/* int	npix		# number of pixels */
/* int	nsinc		# sinc truncation length */
/* int	ntaper		# start of triangular taper */
/* real	staper		# slope of triangular taper */
/* real	mindx		# interpolation minimum */

void c_ii_pcpoly3(float *coeff, int *index, int *len_coeff, float *pcoeff, int *len_pcoeff);
/* real	coeff[ARB]		# 1D array of interpolant coeffcients */
/* int	index			# pointer into coeff array */
/* int	len_coeff		# row length of coeffcients */
/* real	pcoeff[len_pcoeff,ARB]	# polynomial coefficients */
/* int	len_pcoeff		# row length of pcoeff */

void c_ii_pcpoly5(float *coeff, int *index, int *len_coeff, float *pcoeff, int *len_pcoeff);
/* real	coeff[ARB]		# 1D array of interpolant coeffcients */
/* int	index			# pointer into coeff array */
/* int	len_coeff		# row length of coeffcients */
/* real	pcoeff[len_pcoeff,ARB]	# polynomial coefficients */
/* int	len_pcoeff		# row length of pcoeff array */

void c_ii_pcspline3(float *coeff, int *index, int *len_coeff, float *pcoeff, int *len_pcoeff);
/* real	coeff[ARB]		# 1D array of interpolant coeffcients */
/* int	index			# pointer into coeff array */
/* int	len_coeff		# row length of coeffcients */
/* real	pcoeff[len_pcoeff,ARB]	# polynomial coefficients */
/* int	len_pcoeff		# row length of pcoeff */

float c_ii_polterp(float *x, float *y, int *n, float *x0);
/* real	x[ARB],y[ARB]	# x and y array */
/* real	x0		# desired x */
/* int	n		# number of points in x and y = number of */

void c_ii_spline(float *bcoeff, float *diag, int *npts);
/* real	bcoeff[ARB]	# data in and also bspline coefficients out */
/* real	diag[ARB]	# needed for offdiagnol matrix elements */
/* int	npts		# number of data points */

void c_ii_spline2d(float *data, float *coeff, int *nxpix, int *nvectors, int *len_data, int *len_coeff);
/* real	data[len_data,ARB]	# input data array */
/* real	coeff[len_coeff,ARB]	# output array of univariate coefficients in x */
/* int	nxpix			# number of x data points */
/* int	nvectors		# number of univariate splines to calculate */
/* int	len_data		# row dimension of data */
/* int	len_coeff		# row dimension of coeff */

void c_mrider(float *x, float *y, float *datain, int *nxpix, int *nypix, int *len_datain, float *der, int *nxder, int *nyder, int *len_der, int *interp_type);
/* real	x				# x value */
/* real	y				# y value */
/* real	datain[len_datain,ARB]		# data array */
/* int	nxpix				# number of x data points */
/* int	nypix				# number of y data points */
/* int	len_datain			# row length of datain */
/* real	der[len_der, ARB]		# array of derivatives */
/* int	nxder				# number of derivatives in x */
/* int	nyder				# number of derivatives in y */
/* int	len_der				# row length of der, len_der >= nxder */
/* int	interp_type			# interpolant type */

float c_mrieval(float *x, float *y, float *datain, int *nxpix, int *nypix, int *len_datain, int *interp_type);
/* real	x				# x value */
/* real	y				# y value */
/* real	datain[len_datain,ARB]		# data array */
/* int	nxpix				# number of x data points */
/* int	nypix				# number of y data points */
/* int	len_datain			# row length of datain */
/* int	interp_type			# interpolant type */

void c_msider(IRAFPointer *msi, float *x, float *y, float *der, int *nxder, int *nyder, int *len_der);
/* pointer	msi			# pointer to interpolant descriptor structure */
/* real	x			# x value */
/* real	y			# y value */
/* real	der[len_der,ARB]	# derivative array */
/* int	nxder			# number of x derivatives */
/* int	nyder			# number of y derivatives */
/* int	len_der			# row length of der, len_der >= nxder */

float c_msieval(IRAFPointer *msi, float *x, float *y);
/* pointer	msi		# pointer to the interpolant descriptor */
/* real	x		# x data value */
/* real	y		# y data value */

void c_msifit(IRAFPointer *msi, float *datain, int *nxpix, int *nypix, int *len_datain);
/* pointer	msi			# pointer to interpolant descriptor structure */
/* real	datain[len_datain,ARB]	# data array */
/* int	nxpix			# number of points in the x dimension */
/* int	nypix			# number of points in the y dimension */
/* int	len_datain			# row length of datain */

void c_msifree(IRAFPointer *msi);
/* pointer	msi		# pointer to the interpolant descriptor structure */

int c_msigeti(IRAFPointer *msi, int *param);
/* pointer	msi		# interpolant descriptor */
/* int	param		# parameter to be fetched */

void c_msigrid(IRAFPointer *msi, float *x, float *y, float *zfit, int *nx, int *ny, int *len_zfit);
/* pointer	msi 			# pointer to interpolant descriptor structure */
/* real	x[nx]			# array of x values */
/* real	y[ny]			# array of y values */
/* real	zfit[len_zfit,ARB]	# array of fitted values */
/* int	nx			# number of x points */
/* int	ny			# number of y points */
/* int	len_zfit		# row length of zfit */

float c_msigrl(IRAFPointer *msi, float *x, float *y, int *npts);
/* pointer	msi		# pointer to the interpolant descriptor structure */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* int	npts		# number of points which describe the boundary */

void c_ii_find_limits(float *x, float *y, int *npts, float *x1lim, float *x2lim, float *ymin, float *ymax, int *nylmin, int *nylmax);
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* int	npts		# number of data points */
/* real	x1lim[ARB]	# array of x1 limits */
/* real	x2lim[ARB]	# array of x2 limits */
/* real	ymin		# minimum y value for integration */
/* real	ymax		# maximum y value for integration */
/* int	nylmin		# minimum line number for x integration */
/* int	nylmax		# maximum line number for x integration */

int c_ii_pyclip(float *xver, float *yver, float *xintr, float *yintr, int *nver, float *lx, float *ld);
/* real	xver[ARB]		# x vertex coords */
/* real	yver[ARB]		# y vertex coords */
/* real	xintr[ARB]		# x intersection coords */
/* real	yintr[ARB]		# y intersection coords */
/* int	nver			# number of vertices */
/* real	lx, ld 			# equation of image line */

void c_ii_alimrix(float *x, int *npts, float *xmin, int *nxmin, float *xmax, int *nxmax);
/* real	x[npts]		# data */
/* int	npts		# number of data points */
/* real	xmin		# data minimum */
/* int	nxmin		# index of data minimum */
/* real	xmax		# data maximum */
/* int	nxmax		# index of data maximum */

void c_msiinit(IRAFPointer *msi, int *interp_type);
/* pointer	msi		# pointer to the interpolant descriptor structure */
/* int	interp_type	# interpolant type */

void c_msirestore(IRAFPointer *msi, float *interpolant);
/* pointer	msi			# interpolant descriptor */
/* real	interpolant[ARB]	# array containing the interpolant */

void c_msisave(IRAFPointer *msi, float *interpolant);
/* pointer	msi			# interpolant descriptor */
/* real	interpolant[ARB]	# array containing the interpolant */

float c_msisqgrl(IRAFPointer *msi, float *x1, float *x2, float *y1, float *y2);
/* pointer	msi		# pointer to the interpolant descriptor structure */
/* real	x1, x2		# x integration limits */
/* real	y1, y2		# y integration limits */

void c_msivector(IRAFPointer *msi, float *x, float *y, float *zfit, int *npts);
/* pointer	msi		# pointer to the interpolant descriptor structure */
/* real	x[npts]		# array of x values */
/* real	y[npts]		# array of y values */
/* real	zfit[npts]	# array of interpolated values */
/* int	npts		# number of points to be evaluated */

