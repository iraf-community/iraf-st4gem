# include <c_iraf.h>

double c_nlacptsd(IRAFPointer *nl, double *x, double *z, double *w, int *npts, int *nvars);
/* pointer	nl		# pointer to nl fitting structure */
/* double	x[ARB]		# independent variables (npts * nvars) */
/* double	z[ARB]		# function values (npts) */
/* double	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

double c_nlresidd(IRAFPointer *nl, double *x, double *z, double *w, int *npts, int *nvars);
/* pointer	nl		# pointer to nl fitting structure */
/* double	x[ARB]		# independent variables (npts * nvars) */
/* double	z[ARB]		# function values (npts) */
/* double	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

void c_nl_accumd(double *deriv, int *list, double *w, double *dz, int *nfit, double *alpha, double *beta);
/* double	deriv[ARB]	# derivatives */
/* int	list[ARB]	# list of active parameters */
/* double	w		# weight */
/* double	dz		# difference between data and model */
/* int	nfit		# number of fitted parameters */
/* double	alpha[nfit,ARB]	# alpha matrix */
/* double	beta[nfit]	# beta matrix */

float c_nlacptsr(IRAFPointer *nl, float *x, float *z, float *w, int *npts, int *nvars);
/* pointer	nl		# pointer to nl fitting structure */
/* real	x[ARB]		# independent variables (npts * nvars) */
/* real	z[ARB]		# function values (npts) */
/* real	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

float c_nlresidr(IRAFPointer *nl, float *x, float *z, float *w, int *npts, int *nvars);
/* pointer	nl		# pointer to nl fitting structure */
/* real	x[ARB]		# independent variables (npts * nvars) */
/* real	z[ARB]		# function values (npts) */
/* real	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

void c_nl_accumr(float *deriv, int *list, float *w, float *dz, int *nfit, float *alpha, float *beta);
/* real	deriv[ARB]	# derivatives */
/* int	list[ARB]	# list of active parameters */
/* real	w		# weight */
/* real	dz		# difference between data and model */
/* int	nfit		# number of fitted parameters */
/* real	alpha[nfit,ARB]	# alpha matrix */
/* real	beta[nfit]	# beta matrix */

void c_nl_chfacd(double *matrix, int *nbands, int *nrows, double *matfac, int *ier);
/* double   matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* double   matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_nl_chslvd(double *matfac, int *nbands, int *nrows, double *vector, double *coeff);
/* double	matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* double	vector[nrows]			# right side of matrix equation */
/* double	coeff[nrows]			# coefficients */

void c_nl_dampd(double *inmatrix, double *outmatrix, double *constant, int *nbands, int *nrows);
/* double	inmatrix[nbands,ARB]		# input matrix */
/* double	outmatrix[nbands,ARB]		# output matrix */
/* double	constant			# damping constant */
/* int	nbands, nrows			# dimensions of matrix */

void c_nl_chfacr(float *matrix, int *nbands, int *nrows, float *matfac, int *ier);
/* real   matrix[nbands, nrows]	# data matrix */
/* int	nbands			# number of bands */
/* int	nrows			# number of rows */
/* real   matfac[nbands, nrows]	# Cholesky factorization */
/* int	ier			# error code */

void c_nl_chslvr(float *matfac, int *nbands, int *nrows, float *vector, float *coeff);
/* real	matfac[nbands,nrows] 		# Cholesky factorization */
/* int	nbands				# number of bands */
/* int	nrows				# number of rows */
/* real	vector[nrows]			# right side of matrix equation */
/* real	coeff[nrows]			# coefficients */

void c_nl_dampr(float *inmatrix, float *outmatrix, float *constant, int *nbands, int *nrows);
/* real	inmatrix[nbands,ARB]		# input matrix */
/* real	outmatrix[nbands,ARB]		# output matrix */
/* real	constant			# damping constant */
/* int	nbands, nrows			# dimensions of matrix */

void c_nl_dumpd(int *fd, IRAFPointer *nl);
/* int	fd		# file descriptor */
/* pointer	nl		# NLFIT descriptor */

void c_nl_adumpd(int *fd, double *a, int *nrows, int *ncols);
/* int	fd			# file descriptor */
/* double	a[nrows, ncols]		# array */
/* int	nrows, ncols		# dimension */

void c_nl_dumpr(int *fd, IRAFPointer *nl);
/* int	fd		# file descriptor */
/* pointer	nl		# NLFIT descriptor */

void c_nl_adumpr(int *fd, float *a, int *nrows, int *ncols);
/* int	fd			# file descriptor */
/* real	a[nrows, ncols]		# array */
/* int	nrows, ncols		# dimension */

void c_nlerrmsg(int ier, char *errmsg, int maxch);
/* int	ier			# NLFIT error code */
/* char	errmsg[maxch]		# output error message  */
/* int	maxch			# maximum number of chars */

void c_nlerrorsd(IRAFPointer *nl, double *z, double *zfit, double *w, int *npts, double *variance, double *chisqr, double *errors);
/* pointer	nl		# curve descriptor */
/* double	z[ARB]		# data points */
/* double	zfit[ARB]	# fitted data points */
/* double	w[ARB]		# array of weights */
/* int	npts		# number of points */
/* double	variance	# variance of the fit */
/* double	chisqr		# reduced chi-squared of fit (output) */
/* double	errors[ARB]	# errors in coefficients (output) */

void c_nlinvertd(double *alpha, double *chofac, double *cov, double *errors, int *list, int *nfit, double *variance);
/* double	alpha[nfit,ARB]		# data matrix */
/* double	chofac[nfit, ARB]	# cholesky factorization */
/* double	cov[ARB]		# covariance vector */
/* double	errors[ARB]		# computed errors */
/* int	list[ARB]		# list of active parameters */
/* int	nfit			# number of fitted parameters */
/* double	variance		# variance of the fit */

void c_nlerrorsr(IRAFPointer *nl, float *z, float *zfit, float *w, int *npts, float *variance, float *chisqr, float *errors);
/* pointer	nl		# curve descriptor */
/* real	z[ARB]		# data points */
/* real	zfit[ARB]	# fitted data points */
/* real	w[ARB]		# array of weights */
/* int	npts		# number of points */
/* real	variance	# variance of the fit */
/* real	chisqr		# reduced chi-squared of fit (output) */
/* real	errors[ARB]	# errors in coefficients (output) */

void c_nlinvertr(float *alpha, float *chofac, float *cov, float *errors, int *list, int *nfit, float *variance);
/* real	alpha[nfit,ARB]		# data matrix */
/* real	chofac[nfit, ARB]	# cholesky factorization */
/* real	cov[ARB]		# covariance vector */
/* real	errors[ARB]		# computed errors */
/* int	list[ARB]		# list of active parameters */
/* int	nfit			# number of fitted parameters */
/* real	variance		# variance of the fit */

double c_nlevald(IRAFPointer *nl, double *x, int *nvars);
/* pointer	nl		# nlfit descriptor */
/* double	x[ARB]		# x values */
/* int	nvars		# number of variables */

float c_nlevalr(IRAFPointer *nl, float *x, int *nvars);
/* pointer	nl		# nlfit descriptor */
/* real	x[ARB]		# x values */
/* int	nvars		# number of variables */

void c_nlfitd(IRAFPointer *nl, double *x, double *z, double *w, int *npts, int *nvars, int *wtflag, int *stat);
/* pointer	nl		# pointer to nlfit structure */
/* double	x[ARB]		# independent variables (npts * nvars) */
/* double	z[ARB]		# function values (npts) */
/* double	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */
/* int	wtflag		# weighting type */
/* int	stat		# error code */

double c_nlscatterd(IRAFPointer *nl, double *x, double *z, double *w, int *npts, int *nvars);
/* pointer	nl		# Pointer to nl fitting structure */
/* double	x[ARB]		# independent variables (npts * nvars) */
/* double	z[ARB]		# function values (npts) */
/* double	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

void c_nlfitr(IRAFPointer *nl, float *x, float *z, float *w, int *npts, int *nvars, int *wtflag, int *stat);
/* pointer	nl		# pointer to nlfit structure */
/* real	x[ARB]		# independent variables (npts * nvars) */
/* real	z[ARB]		# function values (npts) */
/* real	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */
/* int	wtflag		# weighting type */
/* int	stat		# error code */

float c_nlscatterr(IRAFPointer *nl, float *x, float *z, float *w, int *npts, int *nvars);
/* pointer	nl		# Pointer to nl fitting structure */
/* real	x[ARB]		# independent variables (npts * nvars) */
/* real	z[ARB]		# function values (npts) */
/* real	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

void c_nlfreed(IRAFPointer *nl);
/* pointer	nl		# pointer to non-linear fitting structure */

void c_nlfreer(IRAFPointer *nl);
/* pointer	nl		# pointer to non-linear fitting structure */

void c_nlinitd(IRAFPointer *nl, int *fnc, int *dfnc, double *params, double *dparams, int *nparams, int *plist, int *nfparams, double *tol, int *itmax);
/* pointer	nl		# pointer to nl fitting structure */
/* int	fnc		# fitting function address */
/* int	dfnc		# derivative function address */
/* double	params[ARB]	# initial values for the parameters */
/* double	dparams[ARB]	# initial guess at uncertainties in parameters */
/* int	nparams		# number of parameters */
/* int	plist[ARB]	# list of active parameters */
/* int	nfparams	# number of fitted parameters */
/* double	tol		# fitting tolerance */
/* int	itmax		# maximum number of iterations */

void c_nlinitr(IRAFPointer *nl, int *fnc, int *dfnc, float *params, float *dparams, int *nparams, int *plist, int *nfparams, float *tol, int *itmax);
/* pointer	nl		# pointer to nl fitting structure */
/* int	fnc		# fitting function address */
/* int	dfnc		# derivative function address */
/* real	params[ARB]	# initial values for the parameters */
/* real	dparams[ARB]	# initial guess at uncertainties in parameters */
/* int	nparams		# number of parameters */
/* int	plist[ARB]	# list of active parameters */
/* int	nfparams	# number of fitted parameters */
/* real	tol		# fitting tolerance */
/* int	itmax		# maximum number of iterations */

void c_nliterd(IRAFPointer *nl, double *x, double *z, double *w, int *npts, int *nvars, int *ier);
/* pointer	nl		# pointer to nl fitting structure */
/* double	x[ARB]		# independent variables (npts * nvars) */
/* double	z[ARB]		# function values (npts) */
/* double	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */
/* int	ier		# error code */

void c_nliterr(IRAFPointer *nl, float *x, float *z, float *w, int *npts, int *nvars, int *ier);
/* pointer	nl		# pointer to nl fitting structure */
/* real	x[ARB]		# independent variables (npts * nvars) */
/* real	z[ARB]		# function values (npts) */
/* real	w[ARB]		# weights (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */
/* int	ier		# error code */

void c_nl_list(int *list, int *nlist, int *nfit);
/* int	list[ARB]		# list */
/* int	nlist			# number of elements in the list */
/* int	nfit			# number of active list elments */

void c_nlpgetd(IRAFPointer *nl, double *params, int *nparams);
/* pointer	nl			# pointer to the nlfit structure */
/* double	params[ARB]		# parameter array */
/* int	nparams			# the number of the parameters */

void c_nlpgetr(IRAFPointer *nl, float *params, int *nparams);
/* pointer	nl			# pointer to the nlfit structure */
/* real	params[ARB]		# parameter array */
/* int	nparams			# the number of the parameters */

void c_nlsolved(IRAFPointer *nl, int *ier);
/* pointer	nl	# pointer to the nlfit structure */
/* int	ier	# error code */

void c_nlsolver(IRAFPointer *nl, int *ier);
/* pointer	nl	# pointer to the nlfit structure */
/* int	ier	# error code */

double c_nlstatd(IRAFPointer *nl, int *param);
/* pointer	nl		# pointer to NLFIT structure */
/* int	param		# parameter to be fetched */

int c_nlstati(IRAFPointer *nl, int *param);
/* pointer	nl		# pointer to NLFIT structure */
/* int	param		# parameter to be fetched */

float c_nlstatr(IRAFPointer *nl, int *param);
/* pointer	nl		# pointer to NLFIT structure */
/* int	param		# parameter to be fetched */

void c_nlvectord(IRAFPointer *nl, double *x, double *zfit, int *npts, int *nvars);
/* pointer	nl		# pointer to nl fitting structure */
/* double	x[ARB]		# independent variables (npts * nvars) */
/* double	zfit[ARB]	# function values (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

void c_nlvectorr(IRAFPointer *nl, float *x, float *zfit, int *npts, int *nvars);
/* pointer	nl		# pointer to nl fitting structure */
/* real	x[ARB]		# independent variables (npts * nvars) */
/* real	zfit[ARB]	# function values (npts) */
/* int	npts		# number of points */
/* int	nvars		# number of independent variables */

void c_nlzerod(IRAFPointer *nl);
/* pointer	nl		# pointer to nl fitting structure */

void c_nlzeror(IRAFPointer *nl);
/* pointer	nl		# pointer to nl fitting structure */

