# include <c_iraf.h>

# define IS_UNKNOWN 0
# define IS_IMAGE 1
# define IS_TABLE 2

void c_rdselect(char *file, char *root, char *rowselect, char *colselect, int maxch);
/* char	file[ARB]	# i: filename */
/* char	root[ARB]	# o: filename minus any selectors */
/* char	rowselect[ARB]	# o: row selector */
/* char	colselect[ARB]	# o: column selector */
/* int	maxch		# i: max length of output strings */

void c_tcs_open(IRAFPointer tp, char *columns, IRAFPointer *descrip, int *ndescrip, int maxdescrip);
/* pointer	tp		# i: table descriptor */
/* char	columns[ARB]	# i: list of column names */
/* pointer	descrip[ARB]	# o: list of column array selectors */
/* int	ndescrip	# o: number of column array selectors */
/* int	maxdescrip	# i: length of descrip array */

void c_tcs_close(IRAFPointer *descrip, int ndescrip);
/* pointer	descrip[ARB]	# i: column selectors */
/* int	ndescrip	# i: number of descriptors */

void c_tcs_shape(IRAFPointer descrip, int *length, int *ndim, int maxdimen);
/* pointer	descrip		# i: column selector */
/* int	length[ARB]	# o: dimension lengths */
/* int	ndim		# o: number of dimensions */
/* int	maxdimen	# i: max number of dimensions  */

int c_tcs_totsize(IRAFPointer descrip);
/* pointer	descrip		#i: column selector */

IRAFPointer c_trsopen(IRAFPointer tp, char *expr);
/* pointer	tp		# i: table descriptor */
/* char	expr[ARB]	# i: expression to be parsed */

void c_trsclose(IRAFPointer trs);
/* pointer	trs		# i: Pseudocode structure */

Bool c_trseval(IRAFPointer tp, int irow, IRAFPointer pcode);
/* pointer	tp		# i: table descriptor */
/* int	irow		# i: table row number */
/* pointer	pcode		# i: pseudocode */

IRAFPointer c_rst_copy(IRAFPointer set1);
/* pointer	set1		# i: row set */

IRAFPointer c_rst_create(int loval, int hival);
/* int	loval		# i: low end of range */
/* int	hival		# i: high end of range */

void c_rst_free(IRAFPointer set);
/* pointer	set		# i: row set */

void c_rst_addval(IRAFPointer set, int value);
/* pointer	set		# i: row set */
/* int	value		# i:value to add */

void c_rst_delval(IRAFPointer *set, int value);
/* pointer	set		# u: row set */
/* int	value		# i:value to add */

void c_rst_addtab(IRAFPointer set, int loval, int nval);
/* pointer	set		# i: row set */
/* int	loval		# i: rows are inserted after this row */
/* int	nval		# i: number of rows inserted */

void c_rst_deltab(IRAFPointer *set, int loval, int nval);
/* pointer set		# u: row set */
/* int	loval		# i: first row deleted in underlying table */
/* int	nval		# i: number of rows deleted in underlying table */

IRAFPointer c_rst_and(IRAFPointer set1, IRAFPointer set2);
/* pointer	set1		# i: first row set */
/* pointer	set2		# i: second row set */

IRAFPointer c_rst_or(IRAFPointer set1, IRAFPointer set2);
/* pointer	set1		# i: first row set */
/* pointer	set2		# i: second row set */

IRAFPointer c_rst_not(int nrow, IRAFPointer set1);
/* int	nrow 		# i: largest possible value in set */
/* pointer	set1		# i: set to be negated */

Bool c_rst_inset(IRAFPointer set, int value);
/* pointer	set		# i: row set */
/* int	value		# i: value to be checked */

int c_rst_nelem(IRAFPointer set);
/* pointer	set		# i: row set */

int c_rst_rownum(IRAFPointer set, int index);
/* pointer	set		# i: row set */
/* int	index		# i: index into the set */

void c_rst_show(IRAFPointer set, char *str, int maxch);
/* pointer	set		# i: row set */
/* char	str[ARB]	# o: string representation of set */
/* int	maxch		# i: maximum length of string */

void c_tcs_rdaryb(IRAFPointer tp, IRAFPointer descrip, int irow, int maxbuf, int *nbuf, Bool *buffer);
/* pointer	tp			# i: table descriptor */
/* pointer	descrip			# i: column selector */
/* int	irow			# i: table row number */
/* int	maxbuf			# i: declared length of buffer */
/* int	nbuf			# o: length of output array */
/* bool	buffer[ARB]		# o: array of values */

void c_tcs_rdaryt(IRAFPointer tp, IRAFPointer descrip, int irow, int maxch, int maxbuf, int *nbuf, char **buffer);
/* pointer	tp			# i: table descriptor */
/* pointer	descrip			# i: column selector */
/* int	irow			# i: table row number */
/* int	maxch			# i: max length of string */
/* int	maxbuf			# i: declared length of buffer */
/* int	nbuf			# o: length of output array */
/* char	buffer[maxch,ARB]	# o: array of values */

void c_tcs_rdaryd(IRAFPointer tp, IRAFPointer descrip, int irow, int maxbuf, int *nbuf, double *buffer);
/* pointer	tp			# i: table descriptor */
/* pointer	descrip			# i: column selector */
/* int	irow			# i: table row number */
/* int	maxbuf			# i: declared length of buffer */
/* int	nbuf			# o: length of output array */
/* double	buffer[ARB]		# o: array of values */

void c_tcs_rdaryi(IRAFPointer tp, IRAFPointer descrip, int irow, int maxbuf, int *nbuf, int *buffer);
/* pointer	tp			# i: table descriptor */
/* pointer	descrip			# i: column selector */
/* int	irow			# i: table row number */
/* int	maxbuf			# i: declared length of buffer */
/* int	nbuf			# o: length of output array */
/* int	buffer[ARB]		# o: array of values */

void c_tcs_rdaryr(IRAFPointer tp, IRAFPointer descrip, int irow, int maxbuf, int *nbuf, float *buffer);
/* pointer	tp			# i: table descriptor */
/* pointer	descrip			# i: column selector */
/* int	irow			# i: table row number */
/* int	maxbuf			# i: declared length of buffer */
/* int	nbuf			# o: length of output array */
/* real	buffer[ARB]		# o: array of values */

void c_tcs_rdarys(IRAFPointer tp, IRAFPointer descrip, int irow, int maxbuf, int *nbuf, short *buffer);
/* pointer	tp			# i: table descriptor */
/* pointer	descrip			# i: column selector */
/* int	irow			# i: table row number */
/* int	maxbuf			# i: declared length of buffer */
/* int	nbuf			# o: length of output array */
/* short	buffer[ARB]		# o: array of values */

