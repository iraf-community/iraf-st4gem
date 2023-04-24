# include <c_iraf.h>
IRAFPointer c_immap(char *imspec, int acmode, IRAFPointer hdr_arg);
/* char	imspec[ARB]		#I image specification */
/* int	acmode			#I image access mode */
/* int	hdr_arg			#I length of user fields, or header pointer */

void c_imunmap(IRAFPointer im);
/* pointer	im */

double * c_imgl1d(IRAFPointer im);
/* pointer	im */

int * c_imgl1i(IRAFPointer im);
/* pointer	im */

int * c_imgl1l(IRAFPointer im);
/* pointer	im */

float * c_imgl1r(IRAFPointer im);
/* pointer	im */

short * c_imgl1s(IRAFPointer im);
/* pointer	im */

Complex * c_imgl1x(IRAFPointer im);
/* pointer	im */

double * c_imgl2d(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be read */

int * c_imgl2i(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be read */

int * c_imgl2l(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be read */

float * c_imgl2r(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be read */

short * c_imgl2s(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be read */

Complex * c_imgl2x(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be read */

double * c_imgl3d(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

int * c_imgl3i(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

int * c_imgl3l(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

float * c_imgl3r(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

short * c_imgl3s(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

Complex * c_imgl3x(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

double * c_impl1d(IRAFPointer im);
/* pointer	im			# image header pointer */

int * c_impl1i(IRAFPointer im);
/* pointer	im			# image header pointer */

int * c_impl1l(IRAFPointer im);
/* pointer	im			# image header pointer */

float * c_impl1r(IRAFPointer im);
/* pointer	im			# image header pointer */

short * c_impl1s(IRAFPointer im);
/* pointer	im			# image header pointer */

Complex * c_impl1x(IRAFPointer im);
/* pointer	im			# image header pointer */

double * c_impl2d(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be written */

int * c_impl2i(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be written */

int * c_impl2l(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be written */

float * c_impl2r(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be written */

short * c_impl2s(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be written */

Complex * c_impl2x(IRAFPointer im, int linenum);
/* pointer	im			# image header pointer */
/* int	linenum			# line to be written */

double * c_impl3d(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

int * c_impl3i(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

int * c_impl3l(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

float * c_impl3r(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

short * c_impl3s(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

Complex * c_impl3x(IRAFPointer im, int line, int band);
/* pointer	im		# image header pointer */
/* int	line		# line number within band */
/* int	band		# band number */

int c_imgnld(IRAFPointer imdes, double **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_imgnli(IRAFPointer imdes, int **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_imgnll(IRAFPointer imdes, int **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_imgnlr(IRAFPointer imdes, float **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_imgnls(IRAFPointer imdes, short **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_imgnlx(IRAFPointer imdes, Complex **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_impnld(IRAFPointer imdes, double **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_impnli(IRAFPointer imdes, int **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_impnll(IRAFPointer imdes, int **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_impnlr(IRAFPointer imdes, float **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_impnls(IRAFPointer imdes, short **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

int c_impnlx(IRAFPointer imdes, Complex **lineptr, int *v);
/* pointer	imdes */
/* pointer	lineptr				# on output, points to the pixels */
/* long	v[IM_MAXDIM]			# loop counter */

double * c_imgs1d(IRAFPointer im, int x1, int x2);
/* pointer	im */
/* int	x1, x2 */

int * c_imgs1i(IRAFPointer im, int x1, int x2);
/* pointer	im */
/* int	x1, x2 */

int * c_imgs1l(IRAFPointer im, int x1, int x2);
/* pointer	im */
/* int	x1, x2 */

float * c_imgs1r(IRAFPointer im, int x1, int x2);
/* pointer	im */
/* int	x1, x2 */

short * c_imgs1s(IRAFPointer im, int x1, int x2);
/* pointer	im */
/* int	x1, x2 */

Complex * c_imgs1x(IRAFPointer im, int x1, int x2);
/* pointer	im */
/* int	x1, x2 */

double * c_imgs2d(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

int * c_imgs2i(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

int * c_imgs2l(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

float * c_imgs2r(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

short * c_imgs2s(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

Complex * c_imgs2x(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

double * c_imgs3d(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

int * c_imgs3i(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

int * c_imgs3l(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

float * c_imgs3r(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

short * c_imgs3s(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

Complex * c_imgs3x(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

double * c_imps1d(IRAFPointer im, int x1, int x2);
/* pointer	im		# image header pointer */
/* int	x1		# first column */
/* int	x2		# last column */

int * c_imps1i(IRAFPointer im, int x1, int x2);
/* pointer	im		# image header pointer */
/* int	x1		# first column */
/* int	x2		# last column */

int * c_imps1l(IRAFPointer im, int x1, int x2);
/* pointer	im		# image header pointer */
/* int	x1		# first column */
/* int	x2		# last column */

float * c_imps1r(IRAFPointer im, int x1, int x2);
/* pointer	im		# image header pointer */
/* int	x1		# first column */
/* int	x2		# last column */

short * c_imps1s(IRAFPointer im, int x1, int x2);
/* pointer	im		# image header pointer */
/* int	x1		# first column */
/* int	x2		# last column */

Complex * c_imps1x(IRAFPointer im, int x1, int x2);
/* pointer	im		# image header pointer */
/* int	x1		# first column */
/* int	x2		# last column */

double * c_imps2d(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

int * c_imps2i(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

int * c_imps2l(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

float * c_imps2r(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

short * c_imps2s(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

Complex * c_imps2x(IRAFPointer im, int x1, int x2, int y1, int y2);
/* pointer	im */
/* int	x1, x2, y1, y2 */

double * c_imps3d(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

int * c_imps3i(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

int * c_imps3l(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

float * c_imps3r(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

short * c_imps3s(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

Complex * c_imps3x(IRAFPointer im, int x1, int x2, int y1, int y2, int z1, int z2);
/* pointer	im */
/* int	x1, x2, y1, y2, z1, z2 */

double * c_imggsd(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

int * c_imggsi(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

int * c_imggsl(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

float * c_imggsr(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

short * c_imggss(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

Complex * c_imggsx(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

double * c_impgsd(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

int * c_impgsi(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

int * c_impgsl(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

float * c_impgsr(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

short * c_impgss(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

Complex * c_impgsx(IRAFPointer imdes, int *vs, int *ve, int ndim);
/* pointer	imdes */
/* long	vs[IM_MAXDIM], ve[IM_MAXDIM] */
/* int	ndim */

void c_imflush(IRAFPointer imdes);
/* pointer	imdes */

int c_ximaccess(char *image, int acmode);
/* char	image[ARB]		# image name */
/* int	acmode			# access mode */

void c_imcopy(char *old, char *newname);
/* char	old[ARB]	# old image */
/* char	new[ARB]	# new image */

void c_imdelete(char *image);
/* char	image[ARB] */

void c_imrename(char *old, char *newname);
/* char	old[ARB]		# old image name */
/* char	new[ARB]		# new image name */

void c_imgsection(char *imspec, char *section, int maxch);
/* char	imspec[ARB]		# full image specifcation */
/* char	section[ARB]		# receives image section */
/* int	maxch */

void c_imgimage(char *imspec, char *image, int maxch);
/* char	imspec[ARB]		# full image specification */
/* char	image[ARB]		# receives image name */
/* int	maxch */

void c_imgcluster(char *imspec, char *cluster, int maxch);
/* char	imspec[ARB]		# full image specification */
/* char	cluster[ARB]		# receives root image name */
/* int	maxch */

Bool c_imgetb(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

char c_imgetc(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

double c_imgetd(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

int c_imgeti(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

int c_imgetl(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

float c_imgetr(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

short c_imgets(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */

void c_imgstr(IRAFPointer im, char *key, char *outstr, int maxch);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be returned */
/* char	outstr[ARB]		# output string to receive parameter value */
/* int	maxch */

void c_imputb(IRAFPointer im, char *key, Bool bval);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */
/* bool	bval			# parameter value */

void c_imputd(IRAFPointer im, char *key, double dval);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */
/* double	dval			# double precision value */

void c_imputh(IRAFPointer im, char *key, char *text);
/* pointer	im			#I image descriptor */
/* char	key[ARB]		#I name of the new parameter */
/* char	text[ARB]		#I the history string to be added */

void c_imputi(IRAFPointer im, char *key, int ival);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */
/* int	ival			# parameter value */

void c_imputl(IRAFPointer im, char *key, int lval);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */
/* long	lval			# parameter value */

void c_imputr(IRAFPointer im, char *key, float rval);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */
/* real	rval			# parameter value */

void c_imputs(IRAFPointer im, char *key, short value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */
/* short	value			# parameter value */

void c_impstr(IRAFPointer im, char *key, char *value);
/* pointer	im			#I image descriptor */
/* char	key[ARB]		#I parameter to be set */
/* char	value[ARB]		#I new parameter value */

void c_imaddb(IRAFPointer im, char *key, Bool value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* bool	value			# new or initial value of parameter */

void c_imaddd(IRAFPointer im, char *key, double value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* double	value			# new or initial value of parameter */

void c_imaddi(IRAFPointer im, char *key, int value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* int	value			# new or initial value of parameter */

void c_imaddl(IRAFPointer im, char *key, int value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* long	value			# new or initial value of parameter */

void c_imaddr(IRAFPointer im, char *key, float value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* real	value			# new or initial value of parameter */

void c_imadds(IRAFPointer im, char *key, short value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* short	value			# new or initial value of parameter */

void c_imastr(IRAFPointer im, char *key, char *value);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter or field value */
/* char	value[ARB]		# new or initial value of parameter */

void c_imaddf(IRAFPointer im, char *key, char *datatype);
/* pointer	im			#I image descriptor */
/* char	key[ARB]		#I name of the new parameter */
/* char	datatype[ARB]		#I string permits generalization to domains */

void c_imdelf(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# name of the new parameter */

int c_imaccf(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# name of the new parameter */

int c_imgftype(IRAFPointer im, char *key);
/* pointer	im			# image descriptor */
/* char	key[ARB]		# parameter to be set */

IRAFPointer c_imofnls(IRAFPointer im, char *pattern);
/* pointer	im			# image descriptor */
/* char	template[ARB]		# field name template */

IRAFPointer c_imofnlu(IRAFPointer im, char *pattern);
/* pointer	im			# image descriptor */
/* char	template[ARB]		# field name template */

int c_imgnfn(IRAFPointer fn, char *outstr, int maxch);
/* pointer	fn			# field name list descriptor */
/* char	outstr[ARB]		# output string */
/* int	maxch */

void c_imcfnl(IRAFPointer fn);
/* pointer	fn			# field name list descriptor */

IRAFPointer c_imtopen(char *pattern);
/* char	template[ARB]		# image template */

int c_imtlen(IRAFPointer imt);
/* pointer	imt			# image template descriptor */

void c_imtrew(IRAFPointer imt);
/* pointer	imt			# image template descriptor */

int c_imtgetim(IRAFPointer imt, char *outstr, int maxch);
/* pointer	imt			# image template descriptor */
/* char	outstr[ARB]		# output string */
/* int	maxch			# max chars out */

void c_imtclose(IRAFPointer imt);
/* pointer	imt			# image template descriptor */

IRAFPointer c_imtopenp(char *param);
/* char	param[ARB]		# CL parameter with string value template */

int c_imtrgetim(IRAFPointer imt, int index, char *outstr, int maxch);
/* pointer	imt			# image template descriptor */
/* int	index			# list element to be returned */
/* char	outstr[ARB]		# output string */
/* int	maxch			# max chars out */

int c_imt_mapname(char *fnt, char *outstr, int maxch);
/* char	fnt[ARB]		# FNT string */
/* char	outstr[ARB]		# output string */
/* int	maxch */

int c_imgctime(IRAFPointer imptr);
/* pointer imptr */

void c_imghistory(IRAFPointer imptr, char *str, int maxch);
/* pointer imptr */
/* char 	str[ARB] */
/* int	maxch */

int c_imglimtime(IRAFPointer imptr);
/* pointer imptr */

float c_imgmax(IRAFPointer imptr);
/* pointer imptr */

float c_imgmin(IRAFPointer imptr);
/* pointer imptr */

int c_imgmtime(IRAFPointer imptr);
/* pointer imptr */

int c_imgndim(IRAFPointer imptr);
/* pointer imptr */

int c_imglen(IRAFPointer imptr, int axis);
/* pointer imptr */
/* int axis */

void c_imgpixfile(IRAFPointer imptr, char *str, int maxch);
/* pointer imptr */
/* char 	str[ARB] */
/* int	maxch */

int c_imgtypepix(IRAFPointer imptr);
/* pointer imptr */

void c_imgtitle(IRAFPointer imptr, char *str, int maxch);
/* pointer imptr */
/* char 	str[ARB] */
/* int	maxch */

void c_impctime(IRAFPointer imptr, int val);
/* pointer imptr */
/* long 	val */

void c_imphistory(IRAFPointer imptr, char *str);
/* pointer imptr */
/* char 	str[ARB] */

void c_implimtime(IRAFPointer imptr, int val);
/* pointer imptr */
/* long 	val */

void c_impmax(IRAFPointer imptr, float val);
/* pointer imptr */
/* real 	val */

void c_impmin(IRAFPointer imptr, float val);
/* pointer imptr */
/* real 	val */

void c_impmtime(IRAFPointer imptr, int val);
/* pointer imptr */
/* long 	val */

void c_impndim(IRAFPointer imptr, int val);
/* pointer imptr */
/* int 	val */

void c_implen(IRAFPointer imptr, int axis, int val);
/* pointer imptr */
/* int axis */
/* long 	val */

void c_imppixfile(IRAFPointer imptr, char *str);
/* pointer imptr */
/* char 	str[ARB] */

void c_imptypepix(IRAFPointer imptr, int val);
/* pointer imptr */
/* int 	val */

void c_imptitle(IRAFPointer imptr, char *str);
/* pointer imptr */
/* char 	str[ARB] */

short * c_imguserarea(IRAFPointer imptr);
/* pointer imptr */

