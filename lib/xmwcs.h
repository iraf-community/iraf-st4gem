# include <c_iraf.h>

# define SZ_ATNAME 20

IRAFPointer c_mw_open(IRAFPointer bufptr, int ndim);
/* pointer	bufptr		#I pointer to encoded MWCS, or NULL */
/* int	ndim		#I dimension of system to be created */

IRAFPointer c_mw_openim(IRAFPointer im);
/* pointer	im			#I pointer to image descriptor */

IRAFPointer c_mw_newcopy(IRAFPointer o_mw);
/* pointer	o_mw			#I pointer to old MWCS descriptor */

void c_mw_close(IRAFPointer *mw);
/* pointer	mw		#U pointer to MWCS descriptor */

void c_mw_load(IRAFPointer mw, IRAFPointer bp);
/* pointer	mw			#I pointer to MWCS descriptor */
/* pointer	bp			#I pointer to save buffer, type char */

int c_mw_save(IRAFPointer o_mw, IRAFPointer *bp, int *buflen);
/* pointer	o_mw			#I pointer to MWCS descriptor */
/* pointer	bp			#U pointer to save buffer of type char */
/* int	buflen			#U allocated length of save buffer */

void c_mw_loadim(IRAFPointer mw, IRAFPointer im);
/* pointer	mw			#I pointer to MWCS descriptor */
/* pointer	im			#I pointer to image header */

void c_mw_saveim(IRAFPointer mw, IRAFPointer im);
/* pointer	mw				#I pointer to MWCS descriptor */
/* pointer	im				#I pointer to image descriptor */

IRAFPointer c_mw_sctran(IRAFPointer mw, char *system1, char *system2, int axbits);
/* pointer mw			#I pointer to MWCS descriptor */
/* char	system1[ARB]		#I input coordinate system */
/* char	system2[ARB]		#I output coordinate system */
/* int	axbits			#I bitmap defining axes to be transformed */

int c_mw_gctrand(IRAFPointer a_ct, double *o_ltm, double *o_ltv, int *axtype1, int *axtype2, int maxdim);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* double	o_ltm[ARB]		#O linear tranformation matrix */
/* double	o_ltv[ARB]		#O translation matrix */
/* int	axtype1[ARB]		#O axis types for input system */
/* int	axtype2[ARB]		#O axis types for output system */
/* int	maxdim			#I how much stuff to return */

int c_mw_gctranr(IRAFPointer a_ct, float *o_ltm, float *o_ltv, int *axtype1, int *axtype2, int maxdim);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* real	o_ltm[ARB]		#O linear tranformation matrix */
/* real	o_ltv[ARB]		#O translation matrix */
/* int	axtype1[ARB]		#O axis types for input system */
/* int	axtype2[ARB]		#O axis types for output system */
/* int	maxdim			#I how much stuff to return */

void c_mw_ctfree(IRAFPointer *ct);
/* pointer	ct		#U pointer to CTRAN descriptor */

double c_mw_c1trand(IRAFPointer a_ct, double x);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* double	x			#I coordinates in input system */

float c_mw_c1tranr(IRAFPointer a_ct, float x);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* real	x			#I coordinates in input system */

void c_mw_v1trand(IRAFPointer a_ct, double *x1, double *x2, int npts);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* double	x1[ARB]			#I coordinates in input system */
/* double	x2[ARB]			#O coordinates in output system */
/* int	npts */

void c_mw_v1tranr(IRAFPointer a_ct, float *x1, float *x2, int npts);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* real	x1[ARB]			#I coordinates in input system */
/* real	x2[ARB]			#O coordinates in output system */
/* int	npts */

void c_mw_c2trand(IRAFPointer a_ct, double x1, double y1, double *x2, double *y2);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* double	x1,y1			#I coordinates in input system */
/* double	x2,y2			#O coordinates in output system */

void c_mw_c2tranr(IRAFPointer a_ct, float x1, float y1, float *x2, float *y2);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* real	x1,y1			#I coordinates in input system */
/* real	x2,y2			#O coordinates in output system */

void c_mw_v2trand(IRAFPointer a_ct, double *x1, double *y1, double *x2, double *y2, int npts);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* double	x1[ARB],y1[ARB]		#I coordinates in input system */
/* double	x2[ARB],y2[ARB]		#O coordinates in output system */
/* int	npts */

void c_mw_v2tranr(IRAFPointer a_ct, float *x1, float *y1, float *x2, float *y2, int npts);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* real	x1[ARB],y1[ARB]		#I coordinates in input system */
/* real	x2[ARB],y2[ARB]		#O coordinates in output system */
/* int	npts */

void c_mw_ctrand(IRAFPointer a_ct, double *p1, double *p2, int ndim);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* double	p1[ndim]		#I coordinates of point in input system */
/* double	p2[ndim]		#O coordinates of point in output system */
/* int	ndim			#I dimensionality of point */

void c_mw_ctranr(IRAFPointer a_ct, float *p1, float *p2, int ndim);
/* pointer	a_ct			#I pointer to CTRAN descriptor */
/* real	p1[ndim]		#I coordinates of point in input system */
/* real	p2[ndim]		#O coordinates of point in output system */
/* int	ndim			#I dimensionality of point */

void c_mw_vtrand(IRAFPointer ct, double *v1, double *v2, int ndim, int npts);
/* pointer	ct			#I pointer to CTRAN descriptor */
/* double	v1[ndim,npts]		#I points to be transformed */
/* double	v2[ndim,npts]		#O vector to get the transformed points */
/* int	ndim			#I dimensionality of each point */
/* int	npts			#I number of points */

void c_mw_vtranr(IRAFPointer ct, float *v1, float *v2, int ndim, int npts);
/* pointer	ct			#I pointer to CTRAN descriptor */
/* real	v1[ndim,npts]		#I points to be transformed */
/* real	v2[ndim,npts]		#O vector to get the transformed points */
/* int	ndim			#I dimensionality of each point */
/* int	npts			#I number of points */

void c_mw_sltermd(IRAFPointer mw, double *ltm, double *ltv, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* double	ltm[ndim,ndim]		#I linear transformation matrix */
/* double	ltv[ndim]		#I translation vector */
/* int	ndim			#I dimensionality of system */

void c_mw_sltermr(IRAFPointer mw, float *ltm, float *ltv, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	ltm[ndim,ndim]		#I linear transformation matrix */
/* real	ltv[ndim]		#I translation vector */
/* int	ndim			#I dimensionality of system */

void c_mw_gltermd(IRAFPointer mw, double *ltm, double *ltv, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* double	ltm[ndim,ndim]		#O linear transformation matrix */
/* double	ltv[ndim]		#O translation vector */
/* int	ndim			#I dimensionality of system */

void c_mw_gltermr(IRAFPointer mw, float *ltm, float *ltv, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	ltm[ndim,ndim]		#O linear transformation matrix */
/* real	ltv[ndim]		#O translation vector */
/* int	ndim			#I dimensionality of system */

void c_mw_translated(IRAFPointer mw, double *ltv_1, double *ltm, double *ltv_2, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* double	ltv_1[ndim]		#I input translation vector */
/* double	ltm[ndim,ndim]		#I linear transformation matrix */
/* double	ltv_2[ndim]		#I output translation vector */
/* int	ndim			#I dimensionality of transform */

void c_mw_translater(IRAFPointer mw, float *ltv_1, float *ltm, float *ltv_2, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	ltv_1[ndim]		#I input translation vector */
/* real	ltm[ndim,ndim]		#I linear transformation matrix */
/* real	ltv_2[ndim]		#I output translation vector */
/* int	ndim			#I dimensionality of transform */

void c_mw_rotate(IRAFPointer mw, float theta, float *center, int axbits);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	theta			#I rotation angle, degrees */
/* real	center[ARB]		#I center of rotation */
/* int	axbits			#I bitflags defining axes to be rotated */

void c_mw_scale(IRAFPointer mw, float *scale, int axbits);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	scale[ARB]		#I scale factor for each axis in axbits */
/* int	axbits			#I bitflags defining axes */

void c_mw_shift(IRAFPointer mw, float *shift, int axbits);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	shift[ARB]		#I shift for each axis in axbits */
/* int	axbits			#I bitflags defining axes */

void c_mw_newsystem(IRAFPointer mw, char *system, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* char	system[ARB]		#I system name */
/* int	ndim			#I system dimensionality */

void c_mw_gsystem(IRAFPointer mw, char *outstr, int maxch);
/* pointer	mw			#I pointer to MWCS descriptor */
/* char	outstr[ARB]		#O receives name of world system */
/* int	maxch			#I max chars out */

void c_mw_ssystem(IRAFPointer mw, char *system);
/* pointer	mw			#I pointer to MWCS descriptor */
/* char	system[ARB]		#I system name */

void c_mw_gaxmap(IRAFPointer mw, int *axno, int *axval, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axno[ndim]		#O physical -> logical axis assignments */
/* int	axval[ndim]		#O value of physical axis if axno=0 */
/* int	ndim			#I physical dimension of axis map */

void c_mw_saxmap(IRAFPointer mw, int *axno, int *axval, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axno[ndim]		#I physical -> logical axis assignments */
/* int	axval[ndim]		#I value of physical axis if axno=0 */
/* int	ndim			#I physical dimension of axis map */

void c_mw_gwtermd(IRAFPointer mw, double *r, double *w, double *cd, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* double	r[ndim]			#O physical coordinates of reference point */
/* double	w[ndim]			#O world coordinates of reference point */
/* double	cd[ndim,ndim]		#O CD matrix */
/* int	ndim			#I dimension of Wterm */

void c_mw_gwtermr(IRAFPointer mw, float *r, float *w, float *cd, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	r[ndim]			#O physical coordinates of reference point */
/* real	w[ndim]			#O world coordinates of reference point */
/* real	cd[ndim,ndim]		#O CD matrix */
/* int	ndim			#I dimension of Wterm */

void c_mw_swtermd(IRAFPointer mw, double *r, double *w, double *cd, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* double	r[ndim]			#I physical coordinates of reference point */
/* double	w[ndim]			#I world coordinates of reference point */
/* double	cd[ndim,ndim]		#I CD matrix */
/* int	ndim			#I dimension of Wterm */

void c_mw_swtermr(IRAFPointer mw, float *r, float *w, float *cd, int ndim);
/* pointer	mw			#I pointer to MWCS descriptor */
/* real	r[ndim]			#I physical coordinates of reference point */
/* real	w[ndim]			#I world coordinates of reference point */
/* real	cd[ndim,ndim]		#I CD matrix */
/* int	ndim			#I dimension of Wterm */

void c_mw_swtype(IRAFPointer mw, int *axis, int naxes, char *wtype, char *wattr);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis[naxes]		#I axis number, 1:ndim */
/* int	naxes			#I number of axes in function group */
/* char	wtype[ARB]		#I axis coordinate type */
/* char	wattr[ARB]		#I axis attributes, "attr=value, ..." */

void c_mw_gwsampd(IRAFPointer mw, int axis, double *pv, double *wv, int npts);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis			#I axis which gets the wsamp vector */
/* double	pv[ARB]			#O physical coordinates of points */
/* double	wv[ARB]			#O world coordinates of points */
/* int	npts			#I number of data point in curve */

void c_mw_gwsampr(IRAFPointer mw, int axis, float *pv, float *wv, int npts);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis			#I axis which gets the wsamp vector */
/* real	pv[ARB]			#O physical coordinates of points */
/* real	wv[ARB]			#O world coordinates of points */
/* int	npts			#I number of data point in curve */

void c_mw_swsampd(IRAFPointer mw, int axis, double *pv, double *wv, int npts);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis			#I axis which gets the wsamp vector */
/* double	pv[ARB]			#I physical coordinates of points */
/* double	wv[ARB]			#I world coordinates of points */
/* int	npts			#I number of data point in curve */

void c_mw_swsampr(IRAFPointer mw, int axis, float *pv, float *wv, int npts);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis			#I axis which gets the wsamp vector */
/* real	pv[ARB]			#I physical coordinates of points */
/* real	wv[ARB]			#I world coordinates of points */
/* int	npts			#I number of data point in curve */

void c_mw_gwattrs(IRAFPointer mw, int axis, char *attribute, char *valstr, int maxch);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis			#I axis to which attribute belongs */
/* char	attribute[SZ_ATNAME]	#U attribute name */
/* char	valstr[ARB]		#O attribute value */
/* int	maxch			#I max chars to output value string */

void c_mw_swattrs(IRAFPointer mw, int axis, char *attribute, char *valstr);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	axis			#I axis to which attribute belongs */
/* char	attribute[ARB]		#I attribute name */
/* char	valstr[ARB]		#I attribute value */

void c_mw_invertd(double *o_ltm, double *n_ltm, int ndim);
/* double	o_ltm[ndim,ndim]	#I input matrix */
/* double	n_ltm[ndim,ndim]	#O output (inverted) matrix */
/* int	ndim			#I dimensionality of system */

void c_mw_invertr(float *o_ltm, float *n_ltm, int ndim);
/* real	o_ltm[ndim,ndim]	#I input matrix */
/* real	n_ltm[ndim,ndim]	#O output (inverted) matrix */
/* int	ndim			#I dimensionality of system */

void c_mw_mmuld(double *a, double *b, double *c, int ndim);
/* double	a[ndim,ndim]		#I left input matrix */
/* double	b[ndim,ndim]		#I right input matrix */
/* double	c[ndim,ndim]		#O output matrix */
/* int	ndim			#I dimensionality of system */

void c_mw_mmulr(float *a, float *b, float *c, int ndim);
/* real	a[ndim,ndim]		#I left input matrix */
/* real	b[ndim,ndim]		#I right input matrix */
/* real	c[ndim,ndim]		#O output matrix */
/* int	ndim			#I dimensionality of system */

void c_mw_vmuld(double *a, double *b, double *c, int ndim);
/* double	a[ndim,ndim]		#I input matrix */
/* double	b[ndim]			#I input vector */
/* double	c[ndim]			#O output vector */
/* int	ndim			#I system dimension */

void c_mw_vmulr(float *a, float *b, float *c, int ndim);
/* real	a[ndim,ndim]		#I input matrix */
/* real	b[ndim]			#I input vector */
/* real	c[ndim]			#O output vector */
/* int	ndim			#I system dimension */

void c_mw_seti(IRAFPointer mw, int param, int value);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	param			#I parameter code as defined in <mwset.h> */
/* int	value			#I new value for parameter */

int c_mw_stati(IRAFPointer mw, int param);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	param			#I parameter code as defined in <mwset.h> */

void c_mw_show(IRAFPointer mw, int fd, int what);
/* pointer	mw			#I pointer to MWCS descriptor */
/* int	fd			#I output file */
/* int	what			#I type of output (not used at present) */

