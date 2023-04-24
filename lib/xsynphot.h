# include <c_iraf.h>

void c_compband(char *command, int iw, char *graphtab, char *comptab, int nwave, IRAFPointer *wave, IRAFPointer *band);
/* char	command[ARB]	# i: string containing expression to evaluate */
/* int	iw		# i: position in string where expression starts */
/* char	graphtab[ARB]	# i: graph table name */
/* char	comptab[ARB]	# i: component lookup table name */
/* int	nwave		# i: number of wavelengths */
/* pointer	wave		# u: wavelength array */
/* pointer	band		# o: passband array, calling subroutine must free */

void c_compspec(char *command, int iw, char *graphtab, char *comptab, int nwave, IRAFPointer *wave, IRAFPointer *spec, char *form);
/* char	command[ARB]	# i: string containing expression to evaluate */
/* int	iw		# i: position in string where expression starts */
/* char	graphtab[ARB]	# i: graph table name */
/* char	comptab[ARB]	# i: component lookup table name */
/* int	nwave		# i: number of wavelengths */
/* pointer	wave		# u: wavelength array */
/* pointer	spec		# o: spectrum array, calling array must free */
/* char	form[ARB]	# o: form (units) of spectrum */

void c_inisyntab();

void c_clssyntab();

void c_syncompile(char *expr, int *pcode, int maxcode);
/* char	expr[ARB]	# i: Expression to be parsed */
/* int	pcode[ARB]	# o: Array of pseudocode instructions */
/* int	maxcode		# i: Maximum length of pseudocode array */

void c_syncalc(int *pcode, int maxcode, IRAFPointer getvar, int nwave, float *wave, char *graphtab, char *comptab, float *output, int *units);
/* int	pcode[ARB]	# i: pseudocode used by calculator */
/* int	maxcode		# i: maximum length of pseudocode */
/* pointer	getvar		# i: pointer to subroutine to fetch variable values */
/* int     nwave           # i: length of wavelength and output arrays */
/* real    wave[ARB]       # i: wavelengths at which output is computed */
/* char    graphtab[ARB]   # i: graph table name */
/* char    comptab[ARB]    # i: component lookup table name */
/* real	output[ARB]	# o: result of calculation */
/* int	units		# o: power of FLAM (spectrum = 1, throughput = 0) */

void c_calcrange(int *pcode, int maxcode, char *graphtab, char *comptab, float *minwave, float *maxwave);
/* int	pcode[ARB]	# i: pseudocode used by calculator */
/* int	maxcode		# i: maximum length of pseudocode */
/* char    graphtab[ARB]   # i: graph table name */
/* char    comptab[ARB]    # i: component lookup table name */
/* real	minwave		# o: short end of wavelength range */
/* real	maxwave		# o: long end of wavelength range */

void c_getphotx(char *mode, char *graphtab, char *comptab, char *path, int mxpath, float *phot);
/* char	mode[ARB]	# i: instrument mode */
/* char	graphtab[ARB]	# i: graph table name */
/* char	comptab[ARB]	# i: component lookup table name */
/* char	path[ARB]	# o: list of component throughput tables */
/* int	mxpath		# i: maximum length of path string */
/* real	phot[4]		# o: photometric parameters */

void c_getbandx(char *mode, char *graphtab, char *comptab, Bool logspace, int nwave, float *wave, float *thruput, float *thruerr);
/* char	mode[ARB]	# i: instrument mode */
/* char	graphtab[ARB]	# i: graph table name */
/* char	comptab[ARB]	# i: component lookup table name */
/* bool	logspace	# i: use log spacing for wavelengths? */
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# o: wavelength array */
/* real	thruput[ARB]	# o: grand throughput */
/* real	thruerr[ARB]	# o: throughput error */

void c_evalbandx(char *mode, int nwave, float *wave, char *graphtab, char *comptab, float *thruput, float *thruerr);
/* char	mode[ARB]	# i: instrument mode */
/* int	nwave		# i: length of wavelength and output arrays */
/* real	wave[ARB]	# i: wavelengths at which output is computed */
/* char	graphtab[ARB]	# i: graph table name */
/* char	comptab[ARB]	# i: component lookup table name */
/* real	thruput[ARB]	# o: grand throughput */
/* real	thruerr[ARB]	# o: grand throughput error */

void c_phopar(int nwave, float *wave, float *thruput, float *phot);
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# i: wavelength array */
/* real	thruput[ARB]	# i: throughput array */
/* real	phot[4]		# o: photometric parameters */

int c_anytoang(char *units, float *wave, int nwave);
/* char	units[ARB]	# i: wavelength units */
/* real	wave[ARB]	# u: wavelength array */
/* int	nwave		# i: number of wavelengths */

int c_angtoany(char *units, float *wave, int nwave);
/* char	units[ARB]	# i: wavelength units */
/* real	wave[ARB]	# u: wavelength array */
/* int	nwave		# i: number of wavelengths */

int c_anytophot(char *units, int nwave, float *wave, float *flux);
/* char	units[ARB]	# i: input flux units */
/* int	nwave		# i: length of wavelength and flux arrays */
/* real	wave[ARB]	# i: wavelengths, in angstroms */
/* real	flux[ARB]	# u: flux */

int c_phottoany(char *units, int nwave, float *wave, float *flux);
/* char	units[ARB]	# i: input flux units */
/* int	nwave		# i: length of wavelength and flux arrays */
/* real	wave[ARB]	# i: wavelengths, in angstroms */
/* real	flux[ARB]	# u: flux */

float c_avglam(int nwave, float *wave, float *thruput);
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# i: wavelength array */
/* real	thruput[ARB]	# i: throughput array */

float c_pivlam(int nwave, float *wave, float *thruput);
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# i: wavelength array */
/* real	thruput[ARB]	# i: throughput array */

float c_rmslam(int nwave, float *wave, float *thruput);
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# i: wavelength array */
/* real	thruput[ARB]	# i: throughput array */

float c_fwhmlam(int nwave, float *wave, float *thruput);
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# i: wavelength array */
/* real	thruput[ARB]	# i: throughput array */

float c_funit(float area, int nwave, float *wave, float *thruput);
/* real	area		# i: telescope area */
/* int	nwave		# i: number of wavelengths */
/* real	wave[ARB]	# i: wavelength array */
/* real	thruput[ARB]	# i: throughput array */

void c_getvega(int nwave, float *wave, float *spec);
/* int     nwave           # i: length of wavelength and spectrum arrays */
/* real    wave[ARB]       # i: wavelengths at which spectrum is computed */
/* real	spec[ARB]	# o: spectrum flux */

void c_listpath(char *mode, char *graphtab, char *comptab, char *path, int mxpath);
/* char	mode[ARB]	# i: instrument mode */
/* char	graphtab[ARB]	# i: graph table name */
/* char	comptab[ARB]	# i: component lookup table name */
/* char	path[ARB]	# o: list of component throughput tables */
/* int	mxpath		# i: maximum length of path string */

void c_graffiles(char *mode, char *graphtab, char *comptab, int maxname, int maxthru, int *nthru, char **thruput);
/* char	mode[ARB]			# i: instrument mode */
/* char	graphtab[ARB]			# i: graph table name */
/* char	comptab[ARB]			# i: component lookup table name */
/* int	maxname				# i: max length of thruput table name */
/* int	maxthru				# i: max number of throughput tables */
/* int	nthru				# o: actual number of throughput tables */
/* char	thruput[maxname,maxthru]	# o: throughput table names */

