# include <c_iraf.h>

void c_aabsd(double *a, double *b, int *npix);
/* double	a[ARB], b[ARB] */
/* int	npix, i */

void c_aabsi(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB] */
/* int	npix, i */

void c_aabsl(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB] */
/* int	npix, i */

void c_aabsr(float *a, float *b, int *npix);
/* real	a[ARB], b[ARB] */
/* int	npix, i */

void c_aabss(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB] */
/* int	npix, i */

void c_aabsx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	npix, i */

void c_aaddd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aaddi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aaddkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_aaddki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_aaddkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_aaddkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_aaddks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_aaddkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_aaddl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aaddr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aadds(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aaddx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aandi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aandki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_aandkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_aandks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_aandl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aands(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aavgd(double *a, int *npix, double *mean, double *sigma);
/* double	a[ARB] */
/* int	npix */
/* double	mean, sigma, lcut, hcut */

void c_aavgi(int *a, int *npix, float *mean, float *sigma);
/* int	a[ARB] */
/* int	npix */
/* real	mean, sigma, lcut, hcut */

void c_aavgl(int *a, int *npix, double *mean, double *sigma);
/* long	a[ARB] */
/* int	npix */
/* double	mean, sigma, lcut, hcut */

void c_aavgr(float *a, int *npix, float *mean, float *sigma);
/* real	a[ARB] */
/* int	npix */
/* real	mean, sigma, lcut, hcut */

void c_aavgs(short *a, int *npix, float *mean, float *sigma);
/* short	a[ARB] */
/* int	npix */
/* real	mean, sigma, lcut, hcut */

void c_aavgx(Complex *a, int *npix, float *mean, float *sigma);
/* complex	a[ARB] */
/* int	npix */
/* real	mean, sigma, lcut, hcut */

void c_abavd(double *a, double *b, int *nblocks, int *npix_per_block);
/* double	a[ARB]			# input vector */
/* double	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_abavi(int *a, int *b, int *nblocks, int *npix_per_block);
/* int	a[ARB]			# input vector */
/* int	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_abavl(int *a, int *b, int *nblocks, int *npix_per_block);
/* long	a[ARB]			# input vector */
/* long	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_abavr(float *a, float *b, int *nblocks, int *npix_per_block);
/* real	a[ARB]			# input vector */
/* real	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_abavs(short *a, short *b, int *nblocks, int *npix_per_block);
/* short	a[ARB]			# input vector */
/* short	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_abavx(Complex *a, Complex *b, int *nblocks, int *npix_per_block);
/* complex	a[ARB]			# input vector */
/* complex	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_abeqd(double *a, double *b, int *c, int *npix);
/* double	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abeqi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abeqkd(double *a, double *b, int *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* int	c[ARB] */
/* int	npix */

void c_abeqki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix */

void c_abeqkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* int	c[ARB] */
/* int	npix */

void c_abeqkr(float *a, float *b, int *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* int	c[ARB] */
/* int	npix */

void c_abeqks(short *a, short *b, int *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* int	c[ARB] */
/* int	npix */

void c_abeqkx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* int	c[ARB] */
/* int	npix */

void c_abeql(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abeqr(float *a, float *b, int *c, int *npix);
/* real	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abeqs(short *a, short *b, int *c, int *npix);
/* short	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abeqx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abged(double *a, double *b, int *c, int *npix);
/* double	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgei(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgekd(double *a, double *b, int *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* int	c[ARB] */
/* int	npix */

void c_abgeki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix */

void c_abgekl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* int	c[ARB] */
/* int	npix */

void c_abgekr(float *a, float *b, int *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* int	c[ARB] */
/* int	npix */

void c_abgeks(short *a, short *b, int *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* int	c[ARB] */
/* int	npix */

void c_abgekx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* int	c[ARB] */
/* int	npix */

void c_abgel(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abger(float *a, float *b, int *c, int *npix);
/* real	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abges(short *a, short *b, int *c, int *npix);
/* short	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgex(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgtd(double *a, double *b, int *c, int *npix);
/* double	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgti(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgtkd(double *a, double *b, int *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* int	c[ARB] */
/* int	npix */

void c_abgtki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix */

void c_abgtkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* int	c[ARB] */
/* int	npix */

void c_abgtkr(float *a, float *b, int *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* int	c[ARB] */
/* int	npix */

void c_abgtks(short *a, short *b, int *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* int	c[ARB] */
/* int	npix */

void c_abgtkx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* int	c[ARB] */
/* int	npix */

void c_abgtl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgtr(float *a, float *b, int *c, int *npix);
/* real	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgts(short *a, short *b, int *c, int *npix);
/* short	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abgtx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abled(double *a, double *b, int *c, int *npix);
/* double	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_ablei(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_ablekd(double *a, double *b, int *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* int	c[ARB] */
/* int	npix */

void c_ableki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix */

void c_ablekl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* int	c[ARB] */
/* int	npix */

void c_ablekr(float *a, float *b, int *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* int	c[ARB] */
/* int	npix */

void c_ableks(short *a, short *b, int *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* int	c[ARB] */
/* int	npix */

void c_ablekx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* int	c[ARB] */
/* int	npix */

void c_ablel(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abler(float *a, float *b, int *c, int *npix);
/* real	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_ables(short *a, short *b, int *c, int *npix);
/* short	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_ablex(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abltd(double *a, double *b, int *c, int *npix);
/* double	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_ablti(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abltkd(double *a, double *b, int *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* int	c[ARB] */
/* int	npix */

void c_abltki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix */

void c_abltkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* int	c[ARB] */
/* int	npix */

void c_abltkr(float *a, float *b, int *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* int	c[ARB] */
/* int	npix */

void c_abltks(short *a, short *b, int *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* int	c[ARB] */
/* int	npix */

void c_abltkx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* int	c[ARB] */
/* int	npix */

void c_abltl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abltr(float *a, float *b, int *c, int *npix);
/* real	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_ablts(short *a, short *b, int *c, int *npix);
/* short	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abltx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abned(double *a, double *b, int *c, int *npix);
/* double	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abnei(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abnekd(double *a, double *b, int *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* int	c[ARB] */
/* int	npix */

void c_abneki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix */

void c_abnekl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* int	c[ARB] */
/* int	npix */

void c_abnekr(float *a, float *b, int *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* int	c[ARB] */
/* int	npix */

void c_abneks(short *a, short *b, int *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* int	c[ARB] */
/* int	npix */

void c_abnekx(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* int	c[ARB] */
/* int	npix */

void c_abnel(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abner(float *a, float *b, int *c, int *npix);
/* real	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abnes(short *a, short *b, int *c, int *npix);
/* short	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abnex(Complex *a, Complex *b, int *c, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	c[ARB] */
/* int	npix */

void c_abori(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aborki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_aborkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_aborks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_aborl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_abors(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_absud(double *a, double *b, int *nblocks, int *npix_per_block);
/* double	a[ARB]			# input vector */
/* double	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_absui(int *a, int *b, int *nblocks, int *npix_per_block);
/* int	a[ARB]			# input vector */
/* int	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_absul(int *a, int *b, int *nblocks, int *npix_per_block);
/* long	a[ARB]			# input vector */
/* long	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_absur(float *a, float *b, int *nblocks, int *npix_per_block);
/* real	a[ARB]			# input vector */
/* real	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_absus(short *a, short *b, int *nblocks, int *npix_per_block);
/* short	a[ARB]			# input vector */
/* short	b[nblocks]		# output vector */
/* int	nblocks			# number of blocks (pixels in output vector) */
/* int	npix_per_block		# number of input pixels per block */

void c_achtdd(double *a, double *b, int *npix);
/* double	a[ARB] */
/* double	b[ARB] */
/* int	npix */

void c_achtdi(double *a, int *b, int *npix);
/* double	a[ARB] */
/* int	b[ARB] */
/* int	npix */

void c_achtdl(double *a, int *b, int *npix);
/* double	a[ARB] */
/* long	b[ARB] */
/* int	npix */

void c_achtdr(double *a, float *b, int *npix);
/* double	a[ARB] */
/* real	b[ARB] */
/* int	npix */

void c_achtds(double *a, short *b, int *npix);
/* double	a[ARB] */
/* short	b[ARB] */
/* int	npix */

void c_achtdx(double *a, Complex *b, int *npix);
/* double	a[ARB] */
/* complex	b[ARB] */
/* int	npix */

void c_achtid(int *a, double *b, int *npix);
/* int	a[ARB] */
/* double	b[ARB] */
/* int	npix */

void c_achtii(int *a, int *b, int *npix);
/* int	a[ARB] */
/* int	b[ARB] */
/* int	npix */

void c_achtil(int *a, int *b, int *npix);
/* int	a[ARB] */
/* long	b[ARB] */
/* int	npix */

void c_achtir(int *a, float *b, int *npix);
/* int	a[ARB] */
/* real	b[ARB] */
/* int	npix */

void c_achtis(int *a, short *b, int *npix);
/* int	a[ARB] */
/* short	b[ARB] */
/* int	npix */

void c_achtix(int *a, Complex *b, int *npix);
/* int	a[ARB] */
/* complex	b[ARB] */
/* int	npix */

void c_achtld(int *a, double *b, int *npix);
/* long	a[ARB] */
/* double	b[ARB] */
/* int	npix */

void c_achtli(int *a, int *b, int *npix);
/* long	a[ARB] */
/* int	b[ARB] */
/* int	npix */

void c_achtll(int *a, int *b, int *npix);
/* long	a[ARB] */
/* long	b[ARB] */
/* int	npix */

void c_achtlr(int *a, float *b, int *npix);
/* long	a[ARB] */
/* real	b[ARB] */
/* int	npix */

void c_achtls(int *a, short *b, int *npix);
/* long	a[ARB] */
/* short	b[ARB] */
/* int	npix */

void c_achtlx(int *a, Complex *b, int *npix);
/* long	a[ARB] */
/* complex	b[ARB] */
/* int	npix */

void c_achtrd(float *a, double *b, int *npix);
/* real	a[ARB] */
/* double	b[ARB] */
/* int	npix */

void c_achtri(float *a, int *b, int *npix);
/* real	a[ARB] */
/* int	b[ARB] */
/* int	npix */

void c_achtrl(float *a, int *b, int *npix);
/* real	a[ARB] */
/* long	b[ARB] */
/* int	npix */

void c_achtrr(float *a, float *b, int *npix);
/* real	a[ARB] */
/* real	b[ARB] */
/* int	npix */

void c_achtrs(float *a, short *b, int *npix);
/* real	a[ARB] */
/* short	b[ARB] */
/* int	npix */

void c_achtrx(float *a, Complex *b, int *npix);
/* real	a[ARB] */
/* complex	b[ARB] */
/* int	npix */

void c_achtsd(short *a, double *b, int *npix);
/* short	a[ARB] */
/* double	b[ARB] */
/* int	npix */

void c_achtsi(short *a, int *b, int *npix);
/* short	a[ARB] */
/* int	b[ARB] */
/* int	npix */

void c_achtsl(short *a, int *b, int *npix);
/* short	a[ARB] */
/* long	b[ARB] */
/* int	npix */

void c_achtsr(short *a, float *b, int *npix);
/* short	a[ARB] */
/* real	b[ARB] */
/* int	npix */

void c_achtss(short *a, short *b, int *npix);
/* short	a[ARB] */
/* short	b[ARB] */
/* int	npix */

void c_achtsx(short *a, Complex *b, int *npix);
/* short	a[ARB] */
/* complex	b[ARB] */
/* int	npix */

void c_achtxd(Complex *a, double *b, int *npix);
/* complex	a[ARB] */
/* double	b[ARB] */
/* int	npix */

void c_achtxi(Complex *a, int *b, int *npix);
/* complex	a[ARB] */
/* int	b[ARB] */
/* int	npix */

void c_achtxl(Complex *a, int *b, int *npix);
/* complex	a[ARB] */
/* long	b[ARB] */
/* int	npix */

void c_achtxr(Complex *a, float *b, int *npix);
/* complex	a[ARB] */
/* real	b[ARB] */
/* int	npix */

void c_achtxs(Complex *a, short *b, int *npix);
/* complex	a[ARB] */
/* short	b[ARB] */
/* int	npix */

void c_achtxx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB] */
/* complex	b[ARB] */
/* int	npix */

void c_acjgx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	npix */

void c_aclrd(double *a, int *npix);
/* double	a[ARB] */
/* int	npix, i */

void c_aclri(int *a, int *npix);
/* int	a[ARB] */
/* int	npix, i */

void c_aclrl(int *a, int *npix);
/* long	a[ARB] */
/* int	npix, i */

void c_aclrr(float *a, int *npix);
/* real	a[ARB] */
/* int	npix, i */

void c_aclrs(short *a, int *npix);
/* short	a[ARB] */
/* int	npix, i */

void c_aclrx(Complex *a, int *npix);
/* complex	a[ARB] */
/* int	npix, i */

void c_acnvd(double *in, double *out, int *npix, double *kernel, int *knpix);
/* double	in[npix+knpix-1]	# input vector, including boundary pixels */
/* double	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* double	kernel[knpix]		# convolution kernel */
/* int	knpix			# size of convolution kernel */

void c_acnvi(int *in, int *out, int *npix, int *kernel, int *knpix);
/* int	in[npix+knpix-1]	# input vector, including boundary pixels */
/* int	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* int	kernel[knpix]		# convolution kernel */
/* int	knpix			# size of convolution kernel */

void c_acnvl(int *in, int *out, int *npix, int *kernel, int *knpix);
/* long	in[npix+knpix-1]	# input vector, including boundary pixels */
/* long	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* long	kernel[knpix]		# convolution kernel */
/* int	knpix			# size of convolution kernel */

void c_acnvr(float *in, float *out, int *npix, float *kernel, int *knpix);
/* real	in[npix+knpix-1]	# input vector, including boundary pixels */
/* real	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* real	kernel[knpix]		# convolution kernel */
/* int	knpix			# size of convolution kernel */

void c_acnvrd(double *in, double *out, int *npix, float *kernel, int *knpix);
/* double	in[npix+knpix-1]	# input vector, including boundary pixels */
/* double	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* real	kernel[knpix]		# convolution kernel, always type real */
/* int	knpix			# size of convolution kernel */

void c_acnvri(int *in, int *out, int *npix, float *kernel, int *knpix);
/* int	in[npix+knpix-1]	# input vector, including boundary pixels */
/* int	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* real	kernel[knpix]		# convolution kernel, always type real */
/* int	knpix			# size of convolution kernel */

void c_acnvrl(int *in, int *out, int *npix, float *kernel, int *knpix);
/* long	in[npix+knpix-1]	# input vector, including boundary pixels */
/* long	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* real	kernel[knpix]		# convolution kernel, always type real */
/* int	knpix			# size of convolution kernel */

void c_acnvrr(float *in, float *out, int *npix, float *kernel, int *knpix);
/* real	in[npix+knpix-1]	# input vector, including boundary pixels */
/* real	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* real	kernel[knpix]		# convolution kernel, always type real */
/* int	knpix			# size of convolution kernel */

void c_acnvrs(short *in, short *out, int *npix, float *kernel, int *knpix);
/* short	in[npix+knpix-1]	# input vector, including boundary pixels */
/* short	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* real	kernel[knpix]		# convolution kernel, always type real */
/* int	knpix			# size of convolution kernel */

void c_acnvs(short *in, short *out, int *npix, short *kernel, int *knpix);
/* short	in[npix+knpix-1]	# input vector, including boundary pixels */
/* short	out[ARB]		# output vector */
/* int	npix			# length of output vector */
/* short	kernel[knpix]		# convolution kernel */
/* int	knpix			# size of convolution kernel */

void c_adivd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_adivi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_adivkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_adivki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_adivkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_adivkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_adivks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_adivkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_adivl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_adivr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_adivs(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_adivx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

double c_adotd(double *a, double *b, int *npix);
/* double	a[ARB], b[ARB] */
/* int	npix, i */

float c_adoti(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB] */
/* int	npix, i */

double c_adotl(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB] */
/* int	npix, i */

float c_adotr(float *a, float *b, int *npix);
/* real	a[ARB], b[ARB] */
/* int	npix, i */

float c_adots(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB] */
/* int	npix, i */

float c_adotx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	npix, i */

void c_aexpd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aexpi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aexpkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_aexpki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_aexpkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_aexpkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_aexpks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_aexpkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_aexpl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aexpr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aexps(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aexpx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_afftrr(float *sr, float *si, float *fr, float *fi, int *npix);
/* real	sr[ARB], si[ARB]	# spatial data (input).  SI NOT USED. */
/* real	fr[ARB], fi[ARB]	# real and imag parts of transform (output) */
/* int	npix */

void c_afftrx(float *a, Complex *b, int *npix);
/* real	a[ARB]		# data (input) */
/* complex	b[ARB]		# transform (output).  Dim npix/2+1 */
/* int	npix */

void c_afftxr(float *sr, float *si, float *fr, float *fi, int *npix);
/* real	sr[ARB], si[ARB]	# data, spatial domain (input) */
/* real	fr[ARB], fi[ARB]	# transform, frequency domain (output) */
/* int	npix */

void c_afftxx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB]			# data (input) */
/* complex	b[ARB]			# transform (output) */
/* int	npix */

void c_agltd(double *a, double *b, int *npix, double *low, double *high, double *kmul, double *kadd, int *nrange);
/* double	a[ARB], b[ARB], pixval */
/* int	npix, i */
/* double	low[nrange], high[nrange]	# range limits */
/* double	kmul[nrange], kadd[nrange]	# linear transformation */
/* int	nrange, nr */

void c_aglti(int *a, int *b, int *npix, int *low, int *high, float *kmul, float *kadd, int *nrange);
/* int	a[ARB], b[ARB], pixval */
/* int	npix, i */
/* int	low[nrange], high[nrange]	# range limits */
/* real	kmul[nrange], kadd[nrange] */
/* int	nrange, nr */

void c_agltl(int *a, int *b, int *npix, int *low, int *high, double *kmul, double *kadd, int *nrange);
/* long	a[ARB], b[ARB], pixval */
/* int	npix, i */
/* long	low[nrange], high[nrange]	# range limits */
/* double	kmul[nrange], kadd[nrange]	# linear transformation */
/* int	nrange, nr */

void c_agltr(float *a, float *b, int *npix, float *low, float *high, float *kmul, float *kadd, int *nrange);
/* real	a[ARB], b[ARB], pixval */
/* int	npix, i */
/* real	low[nrange], high[nrange]	# range limits */
/* real	kmul[nrange], kadd[nrange] */
/* int	nrange, nr */

void c_aglts(short *a, short *b, int *npix, short *low, short *high, float *kmul, float *kadd, int *nrange);
/* short	a[ARB], b[ARB], pixval */
/* int	npix, i */
/* short	low[nrange], high[nrange]	# range limits */
/* real	kmul[nrange], kadd[nrange] */
/* int	nrange, nr */

void c_agltx(Complex *a, Complex *b, int *npix, Complex *low, Complex *high, float *kmul, float *kadd, int *nrange);
/* complex	a[ARB], b[ARB], pixval */
/* int	npix, i */
/* complex	low[nrange], high[nrange]	# range limits */
/* real	kmul[nrange], kadd[nrange] */
/* int	nrange, nr */

void c_ahgmd(double *data, int *npix, int *hgm, int *nbins, double *z1, double *z2);
/* double 	data[ARB]		# data vector */
/* int	npix			# number of pixels */
/* int	hgm[ARB]		# output histogram */
/* int	nbins			# number of bins in histogram */
/* double	z1, z2			# greyscale values of first and last bins */

void c_ahgmi(int *data, int *npix, int *hgm, int *nbins, int *z1, int *z2);
/* int 	data[ARB]		# data vector */
/* int	npix			# number of pixels */
/* int	hgm[ARB]		# output histogram */
/* int	nbins			# number of bins in histogram */
/* int	z1, z2			# greyscale values of first and last bins */

void c_ahgml(int *data, int *npix, int *hgm, int *nbins, int *z1, int *z2);
/* long 	data[ARB]		# data vector */
/* int	npix			# number of pixels */
/* int	hgm[ARB]		# output histogram */
/* int	nbins			# number of bins in histogram */
/* long	z1, z2			# greyscale values of first and last bins */

void c_ahgmr(float *data, int *npix, int *hgm, int *nbins, float *z1, float *z2);
/* real 	data[ARB]		# data vector */
/* int	npix			# number of pixels */
/* int	hgm[ARB]		# output histogram */
/* int	nbins			# number of bins in histogram */
/* real	z1, z2			# greyscale values of first and last bins */

void c_ahgms(short *data, int *npix, int *hgm, int *nbins, short *z1, short *z2);
/* short 	data[ARB]		# data vector */
/* int	npix			# number of pixels */
/* int	hgm[ARB]		# output histogram */
/* int	nbins			# number of bins in histogram */
/* short	z1, z2			# greyscale values of first and last bins */

double c_ahivd(double *a, int *npix);
/* double	a[ARB] */
/* int	npix */

int c_ahivi(int *a, int *npix);
/* int	a[ARB] */
/* int	npix */

int c_ahivl(int *a, int *npix);
/* long	a[ARB] */
/* int	npix */

float c_ahivr(float *a, int *npix);
/* real	a[ARB] */
/* int	npix */

short c_ahivs(short *a, int *npix);
/* short	a[ARB] */
/* int	npix */

Complex c_ahivx(Complex *a, int *npix);
/* complex	a[ARB] */
/* int	npix */

void c_aiftrr(float *fr, float *fi, float *sr, float *si, int *npix);
/* real	fr[ARB], fi[ARB]	# real and imag parts of transform (input) */
/* real	sr[ARB], si[ARB]	# spatial data (output).  SI NOT USED. */
/* int	npix */

void c_aiftrx(Complex *a, float *b, int *npix);
/* complex	a[ARB]			# transform, npix/2+1 elements */
/* real	b[ARB]			# output data array */
/* int	npix */

void c_aiftxr(float *fr, float *fi, float *sr, float *si, int *npix);
/* real	fr[ARB], fi[ARB]	# transform, frequency domain (input) */
/* real	sr[ARB], si[ARB]	# data, spatial domain (output) */
/* int	npix */

void c_aiftxx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB]			# transform (input) */
/* complex	b[ARB]			# data (output) */
/* int	npix */

void c_aimgd(Complex *a, double *b, int *npix);
/* complex	a[ARB] */
/* double	b[ARB] */
/* int	npix, i */

void c_aimgi(Complex *a, int *b, int *npix);
/* complex	a[ARB] */
/* int	b[ARB] */
/* int	npix, i */

void c_aimgl(Complex *a, int *b, int *npix);
/* complex	a[ARB] */
/* long	b[ARB] */
/* int	npix, i */

void c_aimgr(Complex *a, float *b, int *npix);
/* complex	a[ARB] */
/* real	b[ARB] */
/* int	npix, i */

void c_aimgs(Complex *a, short *b, int *npix);
/* complex	a[ARB] */
/* short	b[ARB] */
/* int	npix, i */

void c_alimd(double *a, int *npix, double *minval, double *maxval);
/* double	a[ARB], minval, maxval, value */
/* int	npix, i */

void c_alimi(int *a, int *npix, int *minval, int *maxval);
/* int	a[ARB], minval, maxval, value */
/* int	npix, i */

void c_aliml(int *a, int *npix, int *minval, int *maxval);
/* long	a[ARB], minval, maxval, value */
/* int	npix, i */

void c_alimr(float *a, int *npix, float *minval, float *maxval);
/* real	a[ARB], minval, maxval, value */
/* int	npix, i */

void c_alims(short *a, int *npix, short *minval, short *maxval);
/* short	a[ARB], minval, maxval, value */
/* int	npix, i */

void c_alimx(Complex *a, int *npix, Complex *minval, Complex *maxval);
/* complex	a[ARB], minval, maxval, value */
/* int	npix, i */

double c_alovd(double *a, int *npix);
/* double	a[ARB] */
/* int	npix */

int c_alovi(int *a, int *npix);
/* int	a[ARB] */
/* int	npix */

int c_alovl(int *a, int *npix);
/* long	a[ARB] */
/* int	npix */

float c_alovr(float *a, int *npix);
/* real	a[ARB] */
/* int	npix */

short c_alovs(short *a, int *npix);
/* short	a[ARB] */
/* int	npix */

Complex c_alovx(Complex *a, int *npix);
/* complex	a[ARB] */
/* int	npix */

void c_altad(double *a, double *b, int *npix, double *k1, double *k2);
/* double	a[ARB], b[ARB] */
/* double	k1, k2 */
/* int	npix, i */

void c_altai(int *a, int *b, int *npix, float *k1, float *k2);
/* int	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altal(int *a, int *b, int *npix, double *k1, double *k2);
/* long	a[ARB], b[ARB] */
/* double	k1, k2 */
/* int	npix, i */

void c_altar(float *a, float *b, int *npix, float *k1, float *k2);
/* real	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altas(short *a, short *b, int *npix, float *k1, float *k2);
/* short	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altax(Complex *a, Complex *b, int *npix, float *k1, float *k2);
/* complex	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altmd(double *a, double *b, int *npix, double *k1, double *k2);
/* double	a[ARB], b[ARB] */
/* double	k1, k2 */
/* int	npix, i */

void c_altmi(int *a, int *b, int *npix, float *k1, float *k2);
/* int	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altml(int *a, int *b, int *npix, double *k1, double *k2);
/* long	a[ARB], b[ARB] */
/* double	k1, k2 */
/* int	npix, i */

void c_altmr(float *a, float *b, int *npix, float *k1, float *k2);
/* real	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altms(short *a, short *b, int *npix, float *k1, float *k2);
/* short	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altmx(Complex *a, Complex *b, int *npix, float *k1, float *k2);
/* complex	a[ARB], b[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_altrd(double *a, double *b, int *npix, double *k1, double *k2, double *k3);
/* double	a[ARB], b[ARB] */
/* double	k1, k2, k3 */
/* int	npix, i */

void c_altri(int *a, int *b, int *npix, float *k1, float *k2, float *k3);
/* int	a[ARB], b[ARB] */
/* real	k1, k2, k3 */
/* int	npix, i */

void c_altrl(int *a, int *b, int *npix, double *k1, double *k2, double *k3);
/* long	a[ARB], b[ARB] */
/* double	k1, k2, k3 */
/* int	npix, i */

void c_altrr(float *a, float *b, int *npix, float *k1, float *k2, float *k3);
/* real	a[ARB], b[ARB] */
/* real	k1, k2, k3 */
/* int	npix, i */

void c_altrs(short *a, short *b, int *npix, float *k1, float *k2, float *k3);
/* short	a[ARB], b[ARB] */
/* real	k1, k2, k3 */
/* int	npix, i */

void c_altrx(Complex *a, Complex *b, int *npix, float *k1, float *k2, float *k3);
/* complex	a[ARB], b[ARB] */
/* real	k1, k2, k3 */
/* int	npix, i */

void c_aluid(double *a, double *b, float *x, int *npix);
/* double	a[ARB], b[ARB] */
/* real	x[ARB], fraction, tol */
/* int	npix, i, left_pixel */

void c_aluii(int *a, int *b, float *x, int *npix);
/* int	a[ARB], b[ARB] */
/* real	x[ARB], fraction, tol */
/* int	npix, i, left_pixel */

void c_aluil(int *a, int *b, float *x, int *npix);
/* long	a[ARB], b[ARB] */
/* real	x[ARB], fraction, tol */
/* int	npix, i, left_pixel */

void c_aluir(float *a, float *b, float *x, int *npix);
/* real	a[ARB], b[ARB] */
/* real	x[ARB], fraction, tol */
/* int	npix, i, left_pixel */

void c_aluis(short *a, short *b, float *x, int *npix);
/* short	a[ARB], b[ARB] */
/* real	x[ARB], fraction, tol */
/* int	npix, i, left_pixel */

void c_alutd(int *a, double *b, int *npix, double *lut);
/* int	a[ARB]				# input array of indices */
/* double	b[ARB]				# output data array */
/* double	lut[ARB]			# lookup table */
/* int	npix, i */

void c_aluti(int *a, int *b, int *npix, int *lut);
/* int	a[ARB] */
/* int	b[ARB]				# output data array */
/* int	lut[ARB]			# lookup table */
/* int	npix, i */

void c_alutl(int *a, int *b, int *npix, int *lut);
/* long	a[ARB] */
/* long	b[ARB]				# output data array */
/* long	lut[ARB]			# lookup table */
/* int	npix, i */

void c_alutr(int *a, float *b, int *npix, float *lut);
/* int	a[ARB]				# input array of indices */
/* real	b[ARB]				# output data array */
/* real	lut[ARB]			# lookup table */
/* int	npix, i */

void c_aluts(short *a, short *b, int *npix, short *lut);
/* short	a[ARB] */
/* short	b[ARB]				# output data array */
/* short	lut[ARB]			# lookup table */
/* int	npix, i */

void c_amagd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amagi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amagl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amagr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amags(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amagx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amapd(double *a, double *b, int *npix, double *a1, double *a2, double *b1, double *b2);
/* double	a[ARB], b[ARB] */
/* double	a1, a2, b1, b2 */
/* int	npix, i */

void c_amapi(int *a, int *b, int *npix, int *a1, int *a2, int *b1, int *b2);
/* int	a[ARB], b[ARB] */
/* int	a1, a2, b1, b2 */
/* int	npix, i */

void c_amapl(int *a, int *b, int *npix, int *a1, int *a2, int *b1, int *b2);
/* long	a[ARB], b[ARB] */
/* long	a1, a2, b1, b2 */
/* int	npix, i */

void c_amapr(float *a, float *b, int *npix, float *a1, float *a2, float *b1, float *b2);
/* real	a[ARB], b[ARB] */
/* real	a1, a2, b1, b2 */
/* int	npix, i */

void c_amaps(short *a, short *b, int *npix, short *a1, short *a2, short *b1, short *b2);
/* short	a[ARB], b[ARB] */
/* short	a1, a2, b1, b2 */
/* int	npix, i */

void c_amaxd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amaxi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amaxkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_amaxki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_amaxkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_amaxkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_amaxks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_amaxkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_amaxl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amaxr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amaxs(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amaxx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amed3d(double *a, double *b, double *c, double *m, int *npix);
/* double	a[ARB], b[ARB], c[ARB]	# input vectors */
/* double	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed3i(int *a, int *b, int *c, int *m, int *npix);
/* int	a[ARB], b[ARB], c[ARB]	# input vectors */
/* int	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed3l(int *a, int *b, int *c, int *m, int *npix);
/* long	a[ARB], b[ARB], c[ARB]	# input vectors */
/* long	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed3r(float *a, float *b, float *c, float *m, int *npix);
/* real	a[ARB], b[ARB], c[ARB]	# input vectors */
/* real	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed3s(short *a, short *b, short *c, short *m, int *npix);
/* short	a[ARB], b[ARB], c[ARB]	# input vectors */
/* short	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed4d(double *a, double *b, double *c, double *d, double *m, int *npix);
/* double	a[ARB], b[ARB]		# input vectors */
/* double	c[ARB], d[ARB]		# input vectors */
/* double	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed4i(int *a, int *b, int *c, int *d, int *m, int *npix);
/* int	a[ARB], b[ARB]		# input vectors */
/* int	c[ARB], d[ARB]		# input vectors */
/* int	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed4l(int *a, int *b, int *c, int *d, int *m, int *npix);
/* long	a[ARB], b[ARB]		# input vectors */
/* long	c[ARB], d[ARB]		# input vectors */
/* long	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed4r(float *a, float *b, float *c, float *d, float *m, int *npix);
/* real	a[ARB], b[ARB]		# input vectors */
/* real	c[ARB], d[ARB]		# input vectors */
/* real	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed4s(short *a, short *b, short *c, short *d, short *m, int *npix);
/* short	a[ARB], b[ARB]		# input vectors */
/* short	c[ARB], d[ARB]		# input vectors */
/* short	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed5d(double *a, double *b, double *c, double *d, double *e, double *m, int *npix);
/* double	a[ARB], b[ARB]		# input vectors */
/* double	c[ARB], d[ARB], e[ARB]	# input vectors */
/* double	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed5i(int *a, int *b, int *c, int *d, int *e, int *m, int *npix);
/* int	a[ARB], b[ARB]		# input vectors */
/* int	c[ARB], d[ARB], e[ARB]	# input vectors */
/* int	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed5l(int *a, int *b, int *c, int *d, int *e, int *m, int *npix);
/* long	a[ARB], b[ARB]		# input vectors */
/* long	c[ARB], d[ARB], e[ARB]	# input vectors */
/* long	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed5r(float *a, float *b, float *c, float *d, float *e, float *m, int *npix);
/* real	a[ARB], b[ARB]		# input vectors */
/* real	c[ARB], d[ARB], e[ARB]	# input vectors */
/* real	m[ARB]				# output vector (median) */
/* int	npix */

void c_amed5s(short *a, short *b, short *c, short *d, short *e, short *m, int *npix);
/* short	a[ARB], b[ARB]		# input vectors */
/* short	c[ARB], d[ARB], e[ARB]	# input vectors */
/* short	m[ARB]				# output vector (median) */
/* int	npix */

double c_amedd(double *a, int *npix);
/* double	a[ARB] */
/* int	npix */

int c_amedi(int *a, int *npix);
/* int	a[ARB] */
/* int	npix */

int c_amedl(int *a, int *npix);
/* long	a[ARB] */
/* int	npix */

float c_amedr(float *a, int *npix);
/* real	a[ARB] */
/* int	npix */

short c_ameds(short *a, int *npix);
/* short	a[ARB] */
/* int	npix */

Complex c_amedx(Complex *a, int *npix);
/* complex	a[ARB] */
/* int	npix */

void c_amgsd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amgsi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amgsl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amgsr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amgss(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amgsx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amind(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amini(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aminkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_aminki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_aminkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_aminkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_aminks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_aminkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_aminl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aminr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amins(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_aminx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amodd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amodi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amodkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_amodki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_amodkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_amodkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_amodks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_amodl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amodr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amods(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amovd(double *a, double *b, int *npix);
/* double	a[ARB], b[ARB] */
/* int	npix, i, a_first, b_first */

void c_amovi(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB] */
/* int	npix, i, a_first, b_first */

void c_amovkd(double *a, double *b, int *npix);
/* double	a */
/* double	b[ARB] */
/* int	npix, i */

void c_amovki(int *a, int *b, int *npix);
/* int	a */
/* int	b[ARB] */
/* int	npix, i */

void c_amovkl(int *a, int *b, int *npix);
/* long	a */
/* long	b[ARB] */
/* int	npix, i */

void c_amovkr(float *a, float *b, int *npix);
/* real	a */
/* real	b[ARB] */
/* int	npix, i */

void c_amovks(short *a, short *b, int *npix);
/* short	a */
/* short	b[ARB] */
/* int	npix, i */

void c_amovkx(Complex *a, Complex *b, int *npix);
/* complex	a */
/* complex	b[ARB] */
/* int	npix, i */

void c_amovl(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB] */
/* int	npix, i, a_first, b_first */

void c_amovr(float *a, float *b, int *npix);
/* real	a[ARB], b[ARB] */
/* int	npix, i, a_first, b_first */

void c_amovs(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB] */
/* int	npix, i, a_first, b_first */

void c_amovx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	npix, i, a_first, b_first */

void c_amuld(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amuli(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amulkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_amulki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_amulkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_amulkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_amulks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_amulkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_amull(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amulr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amuls(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_amulx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_anegd(double *a, double *b, int *npix);
/* double	a[ARB], b[ARB] */
/* int	npix, i */

void c_anegi(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB] */
/* int	npix, i */

void c_anegl(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB] */
/* int	npix, i */

void c_anegr(float *a, float *b, int *npix);
/* real	a[ARB], b[ARB] */
/* int	npix, i */

void c_anegs(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB] */
/* int	npix, i */

void c_anegx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB] */
/* int	npix, i */

void c_anoti(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB] */
/* int	npix, i */

void c_anotl(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB] */
/* int	npix, i */

void c_anots(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB] */
/* int	npix, i */

void c_apkxd(double *a, double *b, Complex *c, int *npix);
/* double	a[ARB]			# real component */
/* double	b[ARB]			# imaginary component */
/* complex	c[ARB]			# output vector */
/* int	npix, i */

void c_apkxi(int *a, int *b, Complex *c, int *npix);
/* int	a[ARB]			# real component */
/* int	b[ARB]			# imaginary component */
/* complex	c[ARB]			# output vector */
/* int	npix, i */

void c_apkxl(int *a, int *b, Complex *c, int *npix);
/* long	a[ARB]			# real component */
/* long	b[ARB]			# imaginary component */
/* complex	c[ARB]			# output vector */
/* int	npix, i */

void c_apkxr(float *a, float *b, Complex *c, int *npix);
/* real	a[ARB]			# real component */
/* real	b[ARB]			# imaginary component */
/* complex	c[ARB]			# output vector */
/* int	npix, i */

void c_apkxs(short *a, short *b, Complex *c, int *npix);
/* short	a[ARB]			# real component */
/* short	b[ARB]			# imaginary component */
/* complex	c[ARB]			# output vector */
/* int	npix, i */

void c_apkxx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB]			# real component */
/* complex	b[ARB]			# imaginary component */
/* complex	c[ARB]			# output vector */
/* int	npix, i */

double c_apold(double *x, double *coeff, int *ncoeff);
/* double	x		# point at which the polynomial is to be evaluated */
/* double	coeff[ncoeff]	# coefficients of the polynomial, lower orders first */
/* int	ncoeff */

float c_apolr(float *x, float *coeff, int *ncoeff);
/* real	x		# point at which the polynomial is to be evaluated */
/* real	coeff[ncoeff]	# coefficients of the polynomial, lower orders first */
/* int	ncoeff */

void c_apowd(double *a, int *b, double *c, int *npix);
/* double	a[ARB], c[ARB] */
/* int	b[ARB] */
/* int	npix, i */

void c_apowi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], c[ARB] */
/* int	b[ARB] */
/* int	npix, i */

void c_apowkd(double *a, int *b, double *c, int *npix);
/* double	a[ARB], c[ARB] */
/* int	b */
/* int	npix, i */

void c_apowki(int *a, int *b, int *c, int *npix);
/* int	a[ARB], c[ARB] */
/* int	b */
/* int	npix, i */

void c_apowkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], c[ARB] */
/* int	b */
/* int	npix, i */

void c_apowkr(float *a, int *b, float *c, int *npix);
/* real	a[ARB], c[ARB] */
/* int	b */
/* int	npix, i */

void c_apowks(short *a, int *b, short *c, int *npix);
/* short	a[ARB], c[ARB] */
/* int	b */
/* int	npix, i */

void c_apowkx(Complex *a, int *b, Complex *c, int *npix);
/* complex	a[ARB], c[ARB] */
/* int	b */
/* int	npix, i */

void c_apowl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], c[ARB] */
/* int	b[ARB] */
/* int	npix, i */

void c_apowr(float *a, int *b, float *c, int *npix);
/* real	a[ARB], c[ARB] */
/* int	b[ARB] */
/* int	npix, i */

void c_apows(short *a, int *b, short *c, int *npix);
/* short	a[ARB], c[ARB] */
/* int	b[ARB] */
/* int	npix, i */

void c_apowx(Complex *a, int *b, Complex *c, int *npix);
/* complex	a[ARB], c[ARB] */
/* int	b[ARB] */
/* int	npix, i */

int c_aravd(double *a, int *npix, double *mean, double *sigma, double *ksig);
/* double	a[ARB]			# input data array */
/* double	mean, sigma, ksig, deviation, lcut, hcut, lgpx */
/* int	npix, ngpix, old_ngpix, awvgd() */

int c_aravi(int *a, int *npix, float *mean, float *sigma, float *ksig);
/* int	a[ARB]			# input data array */
/* real	mean, sigma, ksig, deviation, lcut, hcut, lgpx */
/* int	npix, ngpix, old_ngpix, awvgi() */

int c_aravl(int *a, int *npix, double *mean, double *sigma, double *ksig);
/* long	a[ARB]			# input data array */
/* double	mean, sigma, ksig, deviation, lcut, hcut, lgpx */
/* int	npix, ngpix, old_ngpix, awvgl() */

int c_aravr(float *a, int *npix, float *mean, float *sigma, float *ksig);
/* real	a[ARB]			# input data array */
/* real	mean, sigma, ksig, deviation, lcut, hcut, lgpx */
/* int	npix, ngpix, old_ngpix, awvgr() */

int c_aravs(short *a, int *npix, float *mean, float *sigma, float *ksig);
/* short	a[ARB]			# input data array */
/* real	mean, sigma, ksig, deviation, lcut, hcut, lgpx */
/* int	npix, ngpix, old_ngpix, awvgs() */

int c_aravx(Complex *a, int *npix, float *mean, float *sigma, float *ksig);
/* complex	a[ARB]			# input data array */
/* real	mean, sigma, ksig, deviation, lcut, hcut, lgpx */
/* int	npix, ngpix, old_ngpix, awvgx() */

void c_arcpd(double *a, double *b, double *c, int *npix);
/* double	a		# constant numerator */
/* double	b[ARB]		# vector denominator */
/* double	c[ARB]		# output vector */
/* int	npix */

void c_arcpi(int *a, int *b, int *c, int *npix);
/* int	a		# constant numerator */
/* int	b[ARB]		# vector denominator */
/* int	c[ARB]		# output vector */
/* int	npix */

void c_arcpl(int *a, int *b, int *c, int *npix);
/* long	a		# constant numerator */
/* long	b[ARB]		# vector denominator */
/* long	c[ARB]		# output vector */
/* int	npix */

void c_arcpr(float *a, float *b, float *c, int *npix);
/* real	a		# constant numerator */
/* real	b[ARB]		# vector denominator */
/* real	c[ARB]		# output vector */
/* int	npix */

void c_arcps(short *a, short *b, short *c, int *npix);
/* short	a		# constant numerator */
/* short	b[ARB]		# vector denominator */
/* short	c[ARB]		# output vector */
/* int	npix */

void c_arcpx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a		# constant numerator */
/* complex	b[ARB]		# vector denominator */
/* complex	c[ARB]		# output vector */
/* int	npix */

void c_argtd(double *a, int *npix, double *ceil, double *newval);
/* double	a[ARB] */
/* int	npix */
/* double	ceil, newval */

void c_argti(int *a, int *npix, int *ceil, int *newval);
/* int	a[ARB] */
/* int	npix */
/* int	ceil, newval */

void c_argtl(int *a, int *npix, int *ceil, int *newval);
/* long	a[ARB] */
/* int	npix */
/* long	ceil, newval */

void c_argtr(float *a, int *npix, float *ceil, float *newval);
/* real	a[ARB] */
/* int	npix */
/* real	ceil, newval */

void c_argts(short *a, int *npix, short *ceil, short *newval);
/* short	a[ARB] */
/* int	npix */
/* short	ceil, newval */

void c_argtx(Complex *a, int *npix, Complex *ceil, Complex *newval);
/* complex	a[ARB] */
/* int	npix */
/* complex	ceil, newval */

void c_arltd(double *a, int *npix, double *floor, double *newval);
/* double	a[ARB] */
/* int	npix */
/* double	floor, newval */

void c_arlti(int *a, int *npix, int *floor, int *newval);
/* int	a[ARB] */
/* int	npix */
/* int	floor, newval */

void c_arltl(int *a, int *npix, int *floor, int *newval);
/* long	a[ARB] */
/* int	npix */
/* long	floor, newval */

void c_arltr(float *a, int *npix, float *floor, float *newval);
/* real	a[ARB] */
/* int	npix */
/* real	floor, newval */

void c_arlts(short *a, int *npix, short *floor, short *newval);
/* short	a[ARB] */
/* int	npix */
/* short	floor, newval */

void c_arltx(Complex *a, int *npix, Complex *floor, Complex *newval);
/* complex	a[ARB] */
/* int	npix */
/* complex	floor, newval */

void c_aseld(double *a, double *b, double *c, int *sel, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b[i] */
/* int	npix */

void c_aseli(int *a, int *b, int *c, int *sel, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b[i] */
/* int	npix */

void c_aselkd(double *a, double *b, double *c, int *sel, int *npix);
/* double	a[ARB], b, c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b */
/* int	npix */

void c_aselki(int *a, int *b, int *c, int *sel, int *npix);
/* int	a[ARB], b, c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b */
/* int	npix */

void c_aselkl(int *a, int *b, int *c, int *sel, int *npix);
/* long	a[ARB], b, c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b */
/* int	npix */

void c_aselkr(float *a, float *b, float *c, int *sel, int *npix);
/* real	a[ARB], b, c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b */
/* int	npix */

void c_aselks(short *a, short *b, short *c, int *sel, int *npix);
/* short	a[ARB], b, c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b */
/* int	npix */

void c_aselkx(Complex *a, Complex *b, Complex *c, int *sel, int *npix);
/* complex	a[ARB], b, c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b */
/* int	npix */

void c_asell(int *a, int *b, int *c, int *sel, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b[i] */
/* int	npix */

void c_aselr(float *a, float *b, float *c, int *sel, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b[i] */
/* int	npix */

void c_asels(short *a, short *b, short *c, int *sel, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b[i] */
/* int	npix */

void c_aselx(Complex *a, Complex *b, Complex *c, int *sel, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	sel[ARB]			# IF sel[i] THEN a[i] ELSE b[i] */
/* int	npix */

double c_asokd(double *a, int *npix, int *ksel);
/* double	a[ARB]			# input array */
/* int	npix			# number of pixels */
/* int	ksel			# element to be selected */

int c_asoki(int *a, int *npix, int *ksel);
/* int	a[ARB]			# input array */
/* int	npix			# number of pixels */
/* int	ksel			# element to be selected */

int c_asokl(int *a, int *npix, int *ksel);
/* long	a[ARB]			# input array */
/* int	npix			# number of pixels */
/* int	ksel			# element to be selected */

float c_asokr(float *a, int *npix, int *ksel);
/* real	a[ARB]			# input array */
/* int	npix			# number of pixels */
/* int	ksel			# element to be selected */

short c_asoks(short *a, int *npix, int *ksel);
/* short	a[ARB]			# input array */
/* int	npix			# number of pixels */
/* int	ksel			# element to be selected */

Complex c_asokx(Complex *a, int *npix, int *ksel);
/* complex	a[ARB]			# input array */
/* int	npix			# number of pixels */
/* int	ksel			# element to be selected */

void c_asrtd(double *a, double *b, int *npix);
/* double	a[ARB], b[ARB]		# input, output arrays */
/* int	npix			# number of pixels */

void c_asrti(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB]		# input, output arrays */
/* int	npix			# number of pixels */

void c_asrtl(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB]		# input, output arrays */
/* int	npix			# number of pixels */

void c_asrtr(float *a, float *b, int *npix);
/* real	a[ARB], b[ARB]		# input, output arrays */
/* int	npix			# number of pixels */

void c_asrts(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB]		# input, output arrays */
/* int	npix			# number of pixels */

void c_asrtx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB]		# input, output arrays */
/* int	npix			# number of pixels */

double c_assqd(double *a, int *npix);
/* double	a[ARB] */
/* int	npix */

float c_assqi(int *a, int *npix);
/* int	a[ARB] */
/* int	npix */

double c_assql(int *a, int *npix);
/* long	a[ARB] */
/* int	npix */

float c_assqr(float *a, int *npix);
/* real	a[ARB] */
/* int	npix */

float c_assqs(short *a, int *npix);
/* short	a[ARB] */
/* int	npix */

Complex c_assqx(Complex *a, int *npix);
/* complex	a[ARB] */
/* int	npix */

void c_asubd(double *a, double *b, double *c, int *npix);
/* double	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_asubi(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_asubkd(double *a, double *b, double *c, int *npix);
/* double	a[ARB] */
/* double	b */
/* double	c[ARB] */
/* int	npix, i */

void c_asubki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_asubkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_asubkr(float *a, float *b, float *c, int *npix);
/* real	a[ARB] */
/* real	b */
/* real	c[ARB] */
/* int	npix, i */

void c_asubks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_asubkx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB] */
/* complex	b */
/* complex	c[ARB] */
/* int	npix, i */

void c_asubl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_asubr(float *a, float *b, float *c, int *npix);
/* real	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_asubs(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_asubx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

double c_asumd(double *a, int *npix);
/* double	a[ARB] */
/* int	npix */

float c_asumi(int *a, int *npix);
/* int	a[ARB] */
/* int	npix */

double c_asuml(int *a, int *npix);
/* long	a[ARB] */
/* int	npix */

float c_asumr(float *a, int *npix);
/* real	a[ARB] */
/* int	npix */

float c_asums(short *a, int *npix);
/* short	a[ARB] */
/* int	npix */

Complex c_asumx(Complex *a, int *npix);
/* complex	a[ARB] */
/* int	npix */

void c_aupxd(Complex *a, double *b, double *c, int *npix);
/* complex	a[ARB]			# input vector */
/* double	b[ARB], c[ARB]	# output vectors */
/* int	npix */

void c_aupxi(Complex *a, int *b, int *c, int *npix);
/* complex	a[ARB]			# input vector */
/* int	b[ARB], c[ARB]	# output vectors */
/* int	npix */

void c_aupxl(Complex *a, int *b, int *c, int *npix);
/* complex	a[ARB]			# input vector */
/* long	b[ARB], c[ARB]	# output vectors */
/* int	npix */

void c_aupxr(Complex *a, float *b, float *c, int *npix);
/* complex	a[ARB]			# input vector */
/* real	b[ARB], c[ARB]	# output vectors */
/* int	npix */

void c_aupxs(Complex *a, short *b, short *c, int *npix);
/* complex	a[ARB]			# input vector */
/* short	b[ARB], c[ARB]	# output vectors */
/* int	npix */

void c_aupxx(Complex *a, Complex *b, Complex *c, int *npix);
/* complex	a[ARB]			# input vector */
/* complex	b[ARB], c[ARB]	# output vectors */
/* int	npix */

Bool c_aveqd(double *a, double *b, int *npix);
/* double	a[ARB], b[ARB]		#I vectors to be compared */
/* int	npix			#I number of pixels to be compared */

Bool c_aveqi(int *a, int *b, int *npix);
/* int	a[ARB], b[ARB]		#I vectors to be compared */
/* int	npix			#I number of pixels to be compared */

Bool c_aveql(int *a, int *b, int *npix);
/* long	a[ARB], b[ARB]		#I vectors to be compared */
/* int	npix			#I number of pixels to be compared */

Bool c_aveqr(float *a, float *b, int *npix);
/* real	a[ARB], b[ARB]		#I vectors to be compared */
/* int	npix			#I number of pixels to be compared */

Bool c_aveqs(short *a, short *b, int *npix);
/* short	a[ARB], b[ARB]		#I vectors to be compared */
/* int	npix			#I number of pixels to be compared */

Bool c_aveqx(Complex *a, Complex *b, int *npix);
/* complex	a[ARB], b[ARB]		#I vectors to be compared */
/* int	npix			#I number of pixels to be compared */

void c_awsud(double *a, double *b, double *c, int *npix, double *k1, double *k2);
/* double	a[ARB], b[ARB], c[ARB] */
/* double	k1, k2 */
/* int	npix, i */

void c_awsui(int *a, int *b, int *c, int *npix, float *k1, float *k2);
/* int	a[ARB], b[ARB], c[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_awsul(int *a, int *b, int *c, int *npix, float *k1, float *k2);
/* long	a[ARB], b[ARB], c[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_awsur(float *a, float *b, float *c, int *npix, float *k1, float *k2);
/* real	a[ARB], b[ARB], c[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_awsus(short *a, short *b, short *c, int *npix, float *k1, float *k2);
/* short	a[ARB], b[ARB], c[ARB] */
/* real	k1, k2 */
/* int	npix, i */

void c_awsux(Complex *a, Complex *b, Complex *c, int *npix, Complex *k1, Complex *k2);
/* complex	a[ARB], b[ARB], c[ARB] */
/* complex	k1, k2 */
/* int	npix, i */

int c_awvgd(double *a, int *npix, double *mean, double *sigma, double *lcut, double *hcut);
/* double	a[ARB] */
/* double	mean, sigma, lcut, hcut */
/* int	npix, i, ngpix */

int c_awvgi(int *a, int *npix, float *mean, float *sigma, float *lcut, float *hcut);
/* int	a[ARB] */
/* real	mean, sigma, lcut, hcut */
/* int	npix, i, ngpix */

int c_awvgl(int *a, int *npix, double *mean, double *sigma, double *lcut, double *hcut);
/* long	a[ARB] */
/* double	mean, sigma, lcut, hcut */
/* int	npix, i, ngpix */

int c_awvgr(float *a, int *npix, float *mean, float *sigma, float *lcut, float *hcut);
/* real	a[ARB] */
/* real	mean, sigma, lcut, hcut */
/* int	npix, i, ngpix */

int c_awvgs(short *a, int *npix, float *mean, float *sigma, float *lcut, float *hcut);
/* short	a[ARB] */
/* real	mean, sigma, lcut, hcut */
/* int	npix, i, ngpix */

int c_awvgx(Complex *a, int *npix, float *mean, float *sigma, float *lcut, float *hcut);
/* complex	a[ARB] */
/* real	mean, sigma, lcut, hcut */
/* int	npix, i, ngpix */

void c_axori(int *a, int *b, int *c, int *npix);
/* int	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_axorki(int *a, int *b, int *c, int *npix);
/* int	a[ARB] */
/* int	b */
/* int	c[ARB] */
/* int	npix, i */

void c_axorkl(int *a, int *b, int *c, int *npix);
/* long	a[ARB] */
/* long	b */
/* long	c[ARB] */
/* int	npix, i */

void c_axorks(short *a, short *b, short *c, int *npix);
/* short	a[ARB] */
/* short	b */
/* short	c[ARB] */
/* int	npix, i */

void c_axorl(int *a, int *b, int *c, int *npix);
/* long	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

void c_axors(short *a, short *b, short *c, int *npix);
/* short	a[ARB], b[ARB], c[ARB] */
/* int	npix, i */

