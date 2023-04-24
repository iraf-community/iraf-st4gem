# include <c_iraf.h>

# define SZ_DATATYPE 79
# define SZ_COMMENT 79
void c_qp_parse(char *qpspec, char *root, int sz_root, char *filter, int sz_filter);
/* char	qpspec[ARB]		#I full event list specification */
/* char	root[sz_root]		#O receives root name */
/* int	sz_root			#I max chars in root name */
/* char	filter[sz_filter]	#O receives filter */
/* int	sz_filter		#I max chars in filter name */

int c_qp_access(char *poefile, int mode);
/* char	poefile[ARB]		#I poefile name */
/* int	mode			#I access mode */

void c_qp_copy(char *o_poefile, char *n_poefile);
/* char	o_poefile[ARB]		#I old poefile name */
/* char	n_poefile[ARB]		#I new poefile name */

void c_qp_rename(char *o_poefile, char *n_poefile);
/* char	o_poefile[ARB]		#I old poefile name */
/* char	n_poefile[ARB]		#I new poefile name */

void c_qp_rebuild(char *poefile);
/* char	poefile[ARB]		#I poefile name */

void c_qp_delete(char *poefile);
/* char	poefile[ARB]		#I poefile name */

IRAFPointer c_qp_open(char *poefile, int mode, IRAFPointer o_qp);
/* char	poefile[ARB]		#I QPOE file to be opened */
/* int	mode			#I file access mode */
/* pointer o_qp			#I reference file, if NEW_COPY */

void c_qp_seti(IRAFPointer qp, int param, int value);
/* pointer	qp			#I QPOE descriptor */
/* int	param			#I parameter to be set */
/* int	value			#I new value for parameter */

int c_qp_stati(IRAFPointer qp, int param);
/* pointer	qp			#I QPOE descriptor */
/* int	param			#I parameter to be queried */

void c_qp_sync(IRAFPointer qp);
/* pointer	qp			#I QPOE descriptor */

void c_qp_close(IRAFPointer qp);
/* pointer	qp			#I QPOE descriptor */

void c_qp_addb(IRAFPointer qp, char *param, Bool value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* bool	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_addc(IRAFPointer qp, char *param, char value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_addd(IRAFPointer qp, char *param, double value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* double	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_addi(IRAFPointer qp, char *param, int value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* int	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_addl(IRAFPointer qp, char *param, int value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* long	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_addr(IRAFPointer qp, char *param, float value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* real	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_adds(IRAFPointer qp, char *param, short value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* short	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_addx(IRAFPointer qp, char *param, Complex value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* complex	value			#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

void c_qp_astr(IRAFPointer qp, char *param, char *value, char *comment);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	value[ARB]		#I parameter value */
/* char	comment[ARB]		#I comment field, if creating parameter */

Bool c_qp_getb(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

Complex c_qp_getx(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

char c_qp_getc(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

double c_qp_getd(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

int c_qp_geti(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

int c_qp_getl(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

float c_qp_getr(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

short c_qp_gets(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

int c_qp_gstr(IRAFPointer qp, char *param, char *outstr, int maxch);
/* pointer qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	outstr[maxch]		#O receives string value */
/* int	maxch			#I max chars out */

void c_qp_putb(IRAFPointer qp, char *param, Bool value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* bool	value			#I scalar parameter value */

void c_qp_putx(IRAFPointer qp, char *param, Complex value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* complex	value			#I scalar parameter value */

void c_qp_putc(IRAFPointer qp, char *param, char value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	value			#I scalar parameter value */

void c_qp_putd(IRAFPointer qp, char *param, double value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* double	value			#I scalar parameter value */

void c_qp_puti(IRAFPointer qp, char *param, int value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* int	value			#I scalar parameter value */

void c_qp_putl(IRAFPointer qp, char *param, int value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* long	value			#I scalar parameter value */

void c_qp_putr(IRAFPointer qp, char *param, float value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* real	value			#I scalar parameter value */

void c_qp_puts(IRAFPointer qp, char *param, short value);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* short	value			#I scalar parameter value */

void c_qp_pstr(IRAFPointer qp, char *param, char *strval);
/* pointer qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	strval[ARB]		#I new string value */

int c_qp_read(IRAFPointer qp, char *param, void *buf, int maxelem, int first, char *datatype);
/* pointer qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	buf[ARB]		#O user data buffer to receive data */
/* int	maxelem			#I max number of data elements to read */
/* int	first			#I first data element to read */
/* char	datatype[ARB]		#I datatype to be returned */

void c_qp_write(IRAFPointer qp, char *param, void *buf, int nelem, int first, char *datatype);
/* pointer qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	buf[ARB]		#I user data buffer containing data */
/* int	nelem			#I number of data elements to write */
/* int	first			#I first data element to write to */
/* char	datatype[ARB]		#I datatype of input data */

int c_qp_popen(IRAFPointer qp, char *param, int mode, int type);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* int	mode			#I file(param) access mode */
/* int	type			#I file type, text or binary */

IRAFPointer c_qp_loadwcs(IRAFPointer qp);
/* pointer	qp				#I QPOE descriptor */

void c_qp_savewcs(IRAFPointer qp, IRAFPointer mw);
/* pointer	qp				#I QPOE descriptor */
/* pointer	mw				#I MWCS descriptor */

int c_qp_accessf(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

void c_qp_deletef(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

void c_qp_renamef(IRAFPointer qp, char *param, char *newname);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	newname[ARB]		#I new parameter name */

void c_qp_copyf(IRAFPointer o_qp, char *o_param, IRAFPointer n_qp, char *n_param);
/* pointer	o_qp			#I QPOE descriptor of old (input) datafile */
/* char	o_param[ARB]		#I input parameter name */
/* pointer	n_qp			#I QPOE descriptor of new (output) datafile */
/* char	n_param[ARB]		#I output parameter name */

void c_qp_addf(IRAFPointer qp, char *param, char *datatype, int maxelem, char *comment, int flags);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	datatype[ARB]		#I parameter data type */
/* int	maxelem			#I allocated length of parameter */
/* char	comment[ARB]		#I comment describing parameter */
/* int	flags			#I parameter flags */

int c_qp_queryf(IRAFPointer qp, char *param, char *datatype, int *maxelem, char *comment, int *flags);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */
/* char	datatype[SZ_DATATYPE]	#O parameter data type */
/* int	maxelem			#O allocated length of parameter */
/* char	comment[SZ_COMMENT]	#O comment describing parameter */
/* int	flags			#O parameter flag word */

int c_qp_lenf(IRAFPointer qp, char *param);
/* pointer	qp			#I QPOE descriptor */
/* char	param[ARB]		#I parameter name */

int c_qp_expandtext(IRAFPointer qp, char *s1, char *s2, int maxch);
/* pointer	qp			#I QPOE descriptor */
/* char	s1[ARB]			#I input string containing macros */
/* char	s2[maxch]		#O output string buffer */
/* int	maxch			#I max chars out */

IRAFPointer c_qp_ofnls(IRAFPointer qp, char *fieldname);
/* pointer	qp			#I QPOE descriptor */
/* char	template[ARB]		#I field name template */

IRAFPointer c_qp_ofnlu(IRAFPointer qp, char *fieldname);
/* pointer	qp			#I QPOE descriptor */
/* char	template[ARB]		#I field name template */

IRAFPointer c_qp_ofnl(IRAFPointer qp, char *fieldname, Bool sort);
/* pointer	qp			#I QPOE descriptor */
/* char	template[ARB]		#I field name template */
/* bool	sort			#I sort list of matched names? */

int c_qp_gnfn(IRAFPointer fl, char *outstr, int maxch);
/* pointer	fl			#I list descriptor */
/* char	outstr[maxch]		#O output string */
/* int	maxch			#I max chars out */

int c_qp_lenfnl(IRAFPointer fl);
/* pointer	fl			#I list descriptor */

void c_qp_seekfnl(IRAFPointer fl, int pos);
/* pointer	fl			#I list descriptor */
/* int	pos			#I desired list element, BOF, EOF */

void c_qp_cfnl(IRAFPointer fl);
/* pointer	fl			#I list descriptor */

IRAFPointer c_qpio_open(IRAFPointer qp, char *paramex, int mode);
/* pointer	qp			#I QPOE descriptor */
/* char	paramex[ARB]		#I event-list parameter plus expression list */
/* int	mode			#I access mode */

void c_qpio_close(IRAFPointer io);
/* pointer	io			#I QPIO descriptor */

void c_qpio_setrange(IRAFPointer io, int *vs, int *ve, int ndim);
/* pointer	io			#I QPIO descriptor */
/* int	vs[ARB]			#I start vector (lower left corner) */
/* int	ve[ARB]			#I end vector (upper right corner) */
/* int	ndim			#I vector length (ndim=2 at present) */

int c_qpio_getrange(IRAFPointer io, int *vs, int *ve, int maxdim);
/* pointer	io			#I QPIO descriptor */
/* int	vs[ARB]			#O start vector (lower left corner) */
/* int	ve[ARB]			#O end vector (upper right corner) */
/* int	maxdim			#I vector length (ndim=2 at present) */

void c_qpio_setfilter(IRAFPointer io, char *expr);
/* pointer	io			#I QPIO descriptor */
/* char	expr[ARB]		#I option setting expression */

int c_qpio_getfilter(IRAFPointer io, char *outstr, int maxch);
/* pointer	io			#I QPIO descriptor */
/* char	outstr[maxch]		#O where to put the filter text */
/* int	maxch			#I max chars out */

void c_qpio_seti(IRAFPointer io, int param, int value);
/* pointer	io			#I QPIO descriptor */
/* int	param			#I parameter code */
/* int	value			#I new parameter value */

int c_qpio_stati(IRAFPointer io, int param);
/* pointer	io			#I QPIO descriptor */
/* int	param			#I parameter code */

IRAFPointer c_qpio_loadwcs(IRAFPointer io);
/* pointer	io				#I QPIO descriptor */

void c_qpio_mkindex(IRAFPointer io, char *key);
/* pointer	io			#I QPIO descriptor */
/* char	key[ARB]		#I list of key fields */

void c_qpio_putevents(IRAFPointer io, IRAFPointer *i_ev, int nevents);
/* pointer	io			#I QPIO descriptor */
/* pointer	i_ev[ARB]		#I array of event pointers */
/* int	nevents			#I number of events */

int c_qpio_getevents(IRAFPointer io, IRAFPointer *o_ev, int *maskval, int maxev, int *o_nev);
/* pointer	io			#I QPIO descriptor */
/* pointer	o_ev[maxev]		#O receives the event struct pointers */
/* int	maskval			#O receives the mask value of the events */
/* int	maxev			#I max events out */
/* int	o_nev			#O same as function value (nev_out|EOF) */

void c_qpio_put(IRAFPointer io, char **buffer, int numrows);
/* pointer	io			#I QPIO descriptor */
/* char         **buffer		#I buffer of char **data */
/* int	nevents			#I number of events */

int c_qpio_get(IRAFPointer io, int buffer_size, char **IrafRowPtrs, 
               int maskval, int maxev, int o_nev); 
/* pointer	io			#I QPIO descriptor */
/* int  buffer_size             #I number of IRAF pointers to IRAF memory */
/* char **IrafRowPtrs           #I/O IRAF ptrs to memory converted to char* */
/* int	maskval			#O receives the mask value of the events */
/* int	maxev			#I max events out */
/* int	o_nev			#O same as function value (nev_out|EOF) */

int c_qpio_readpixi(IRAFPointer io, int *obuf, int *vs, int *ve, int ndim, int xblock, int yblock);
/* pointer	io			#I QPIO descriptor */
/* int	obuf[ARB]		#O output pixel buffer */
/* int	vs[ndim], ve[ndim]	#I vectors defining region to be extracted */
/* int	ndim			#I should be 2 for QPOE */
/* int	xblock, yblock		#I blocking factors */

int c_qpio_readpixs(IRAFPointer io, short *obuf, int *vs, int *ve, int ndim, int xblock, int yblock);
/* pointer	io			#I QPIO descriptor */
/* short	obuf[ARB]		#O output pixel buffer */
/* int	vs[ndim], ve[ndim]	#I vectors defining region to be extracted */
/* int	ndim			#I should be 2 for QPOE */
/* int	xblock, yblock		#I blocking factors */

IRAFPointer c_qpex_open(IRAFPointer qp, char *expr);
/* pointer	qp			#I QPOE descriptor */
/* char	expr[ARB]		#I selection expression (filter)  */

int c_qpex_modfilter(IRAFPointer ex, char *exprlist);
/* pointer ex                      #I qpex descriptor */
/* char    exprlist[ARB]           #I list of attribute=expr expressions */

int c_qpex_getfilter(IRAFPointer ex, char *outstr, int maxch);
/* pointer	ex			#I QPEX descriptor */
/* char	outstr[maxch]		#O receives the filter string */
/* int	maxch			#I max chars out */

int c_qpex_getattribute(IRAFPointer ex, char *attribute, char *outstr, int maxch);
/* pointer	ex			#I QPEX descriptor */
/* char	attribute[ARB]		#I attribute name */
/* char	outstr[maxch]		#O receives the filter string */
/* int	maxch			#I max chars out */

int c_qpex_attrld(IRAFPointer ex, char *attribute, IRAFPointer *xs, IRAFPointer *xe, int *xlen);
/* pointer	ex			#I QPEX descriptor */
/* char	attribute[ARB]		#I attribute name */
/* pointer	xs			#U pointer to array of start values */
/* pointer	xe			#U pointer to array of end values */
/* int	xlen			#U length of xs/xe arrays */

int c_qpex_attrli(IRAFPointer ex, char *attribute, IRAFPointer *xs, IRAFPointer *xe, int *xlen);
/* pointer	ex			#I QPEX descriptor */
/* char	attribute[ARB]		#I attribute name */
/* pointer	xs			#U pointer to array of start values */
/* pointer	xe			#U pointer to array of end values */
/* int	xlen			#U length of xs/xe arrays */

int c_qpex_attrlr(IRAFPointer ex, char *attribute, IRAFPointer *xs, IRAFPointer *xe, int *xlen);
/* pointer	ex			#I QPEX descriptor */
/* char	attribute[ARB]		#I attribute name */
/* pointer	xs			#U pointer to array of start values */
/* pointer	xe			#U pointer to array of end values */
/* int	xlen			#U length of xs/xe arrays */

int c_qpex_evaluate(IRAFPointer ex, IRAFPointer *i_ev, IRAFPointer *o_ev, int nev);
/* pointer ex			#I QPEX descriptor (expression) */
/* pointer i_ev[nev]		#I array of pointers to event structs */
/* pointer	o_ev[nev]		#O receives the pointers of the passed events */
/* int	nev			#I number of input events */

void c_qpex_close(IRAFPointer ex);
/* pointer	ex			#I QPEX descriptor */

