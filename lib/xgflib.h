# include <c_iraf.h>
# define SZ_KEYWORD  8
# define SZ_SHORTSTR 32
# define TYPE_LIST   "|stf|fxf|oif|plf|qpf|"
# define UNKNOWN_FMT 0
# define GEIS_FMT    1
# define FITS_FMT    2
# define IRAF_FMT    3
# define PIXLIST_FMT 4
# define QPOE_FMT    5
# define KEYSTRING   "|gcount|psize|pcount|group|bitpix|szgroup|"
# define GCOUNT_KEY  1
# define PSIZE_KEY   2
# define PCOUNT_KEY  3
# define GROUP_KEY   4
# define BITPIX_KEY  5
# define SZGROUP_KEY 6

void c_gf_addpar(IRAFPointer im, char *pname, int dtype, int plen, char *pval, char *pcomm);
/* pointer	im			# i: image descriptor */
/* char	pname[ARB]		# i: parameter name */
/* int	dtype			# i: SPP datatype of parameter */
/* int	plen			# i: length (> 1 if array) */
/* char	pval[ARB]		# i: string encoded initial parameter value */
/* char	pcomm[ARB]		# i: string comment to the new parameter */

void c_gf_delpar(IRAFPointer im, char *pname);
/* pointer	im		# i: image descriptor */
/* char	pname[ARB]	# i: parameter name */

int c_gf_gcount(IRAFPointer im);
/* pointer	im		# i: image descriptor */

Bool c_gf_geis(IRAFPointer im);
/* pointer	im		# i: Image descriptor */

int c_gf_gfind(IRAFPointer im, char *keyword);
/* pointer im              # i: Image descriptor */
/* char    keyword[ARB]    # i: Group parameter keyword name */

int c_gf_gstfval(IRAFPointer im, char *keyword);
/* pointer	im		# i: image descriptor */
/* char	keyword[ARB]	# i: keyword to retrieve */

IRAFPointer c_gf_map(char *image, int acmode, IRAFPointer oldim);
/* char	image[ARB]	# i: image name */
/* int	acmode		# i: access mode */
/* pointer	oldim		# i: old image descriptor */

void c_gf_opengr(IRAFPointer *im, int gn, float *datamin, float *datamax, IRAFPointer imt);
/* pointer	im		# u: image descriptor */
/* int	gn		# i: group number to skip to */
/* real	datamin		# u: image minimun value */
/* real	datamax		# u: image maximum value */
/* pointer imt		# i: image template descriptor (NEW_COPY only) */

void c_gf_pstfval(IRAFPointer im, char *keyword, int value);
/* pointer	im		# i: image descriptor */
/* char	keyword[ARB]	# i: keyword name */
/* int	value		# i: keyword value */

