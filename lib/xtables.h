# include <c_iraf.h>

/* These constants do NOT include space for the NULL character */
# define SZ_KEYWORD 8
# define SZ_COLNAME 79
# define SZ_COLUNITS 79
# define SZ_COLFMT 79
# define SZ_PARREC 80

/* Indefinite valued numbers. Taken from IRAF's spp.h */
# define IRAF_INDEFS (-32767)
# define IRAF_INDEFL (0x80000001)
# define IRAF_INDEFI IRAF_INDEFL
# define IRAF_INDEFR 1.6e38
# if defined(__VMS)
# define IRAF_INDEFD 1.6e38
# else
# define IRAF_INDEFD 1.6e308
# endif

/* These may be set by tbpset and/or read by tbpsta: */
/* (R) row length to allocate (SZ_REAL) */
# define TBL_ROWLEN  1 
/* (R) increase row length (SZ_REAL) */
# define TBL_INCR_ROWLEN  2 
/* (C) number of rows to allocate */
# define TBL_ALLROWS  3 
/* (C) increase alloc num of rows */
# define TBL_INCR_ALLROWS 4 
/* which type of table? */
# define TBL_WHTYPE  5 
/* maximum number of user param */
# define TBL_MAXPAR  6 
/* maximum number of columns */
# define TBL_MAXCOLS  7 
/* (R) row length to allocate (SZ_CHAR) */
# define TBL_ROWLEN_CHAR  8 
/* machine independent table format */
# define TBL_TYPE_MI  10 
/* row-ordered table format */
# define TBL_TYPE_S_ROW  11 
/* column-ordered table format */
# define TBL_TYPE_S_COL  12 
/* text file */
# define TBL_TYPE_TEXT  13 
/* FITS table */
# define TBL_TYPE_FITS  14 
/* don't use this yet! */
# define TBL_TYPE_CDF  15 
/* These may be read by tbpsta but may not be set: */
/* number of rows written to */
# define TBL_NROWS  21 
/* number of columns # defined */
# define TBL_NCOLS  22 
/* (R) amount of row len used (SZ_REAL) */
# define TBL_ROWLEN_USED  23 
/* number of user parameters */
# define TBL_NPAR  24 
/* (R) amount of row len used (SZ_CHAR) */
# define TBL_ROWLEN_CHAR_USED 25 
/* version that created the table */
# define TBL_VERSION  26 
/* These have to do with the file size and file I/O buffer size. */
/* set RANDOM or SEQUENTIAL */
# define TBL_ADVICE  31 
/* get buffer size in char */
# define TBL_BUFSIZE  32 
/* get size of table data in char */
# define TBL_DATA_SIZE  34 
/* These are for information about a column. */
/* column name */
# define TBL_COL_NAME  41 
/* units for column */
# define TBL_COL_UNITS  42 
/* print format for displaying values */
# define TBL_COL_FMT  43 
/* data type (-n for char string) */
# define TBL_COL_DATATYPE 44 
/* column number */
# define TBL_COL_NUMBER  45 
/* length for printing using print fmt */
# define TBL_COL_FMTLEN  46 
/* number of elements if it's an array */
# define TBL_COL_LENDATA  47 
/* dimension of array */
# define TBL_COL_DIMENSION 48

IRAFPointer c_tbtopn(char *tablename, int iomode, IRAFPointer oldname);
/* char	tablename[ARB]	# i: the name of the table */
/* int	iomode		# i: I/O mode */
/* pointer template	# i: pointer to template table, or zero */

void c_tbtcre(IRAFPointer tp);
/* pointer tp			# Pointer to table descriptor */

void c_tbtclo(IRAFPointer tp);
/* pointer tp			# i: pointer to descriptor of table to be closed */

void c_tbcdef(IRAFPointer tp, IRAFPointer *colptr, char **colname, char **colunits, char **colfmt, int *datatype, int *lenarray, int numcols);
/* pointer tp				# i: pointer to table descriptor */
/* pointer colptr[numcols]			# o: pointers to new columns */
/* char	colname[SZ_COLNAME,numcols]	# i: names of columns */
/* char	colunits[SZ_COLUNITS,numcols]	# i: units for columns */
/* char	colfmt[SZ_COLFMT,numcols]	# i: print formats for columns */
/* int	datatype[numcols]		# i: data types of columns */
/* int	lenarray[numcols]		# i: number of elements for each column */
/* int	numcols				# i: number of columns to be defined */

void c_tbcdef1(IRAFPointer tp, IRAFPointer *colptr, char *colname, char *colunits, char *colfmt, int datatype, int nelem);
/* pointer tp			# i: pointer to table descriptor */
/* pointer colptr			# o: pointer to new column */
/* char	colname[SZ_COLNAME]	# i: name of column */
/* char	colunits[SZ_COLUNITS]	# i: units for column */
/* char	colfmt[SZ_COLFMT]	# i: print format for column */
/* int	datatype		# i: data types of column */
/* int	nelem			# i: number of elements for column */

void c_tbcfnd(IRAFPointer tp, char **colname, IRAFPointer *colptr, int numcols);
/* pointer tp				# i: pointer to table descriptor */
/* char	colname[SZ_COLNAME,numcols]	# i: array of column names */
/* pointer colptr[numcols]			# o: array of ptrs to column descriptors */
/* int	numcols				# i: length of arrays colname & colptr */

void c_tbcfnd1(IRAFPointer tp, char *colname, IRAFPointer *colptr);
/* pointer tp			# i: pointer to table descriptor */
/* char	colname[SZ_COLNAME]	# i: name of column */
/* pointer colptr			# o: pointer to column, or NULL */

void c_tbcinf(IRAFPointer colptr, int *colnum, char *colname, char *colunits, char *colfmt, int *datatype, int *lendata, int *lenfmt);
/* pointer colptr			# i: Pointer to a column descriptor */
/* int	colnum			# o: Column number */
/* char	colname[ARB]		# o: Column name */
/* char	colunits[ARB]		# o: Units for column */
/* char	colfmt[ARB]		# o: Print format for display of column */
/* int	datatype		# o: Data type of column (SPP type or -n) */
/* int	lendata			# o: Number of elements (=1) */
/* int	lenfmt			# o: Bytes for print format */

int c_tbcigi(IRAFPointer cptr, int get_what);
/* pointer cptr			# i: pointer to column descriptor */
/* int	get_what		# i: indicates what column info to get */

void c_tbciga(IRAFPointer tp, IRAFPointer cp, int *ndim, int *axlen, int maxdim);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to column descriptor */
/* int	ndim		# o: dimension of array */
/* int	axlen[maxdim]	# o: length of each axis */
/* int	maxdim		# i: size of axlen array */

void c_tbcigt(IRAFPointer cptr, int get_what, char *outstr, int maxch);
/* pointer cptr			# i: pointer to column descriptor */
/* int	get_what		# i: indicates what string to get */
/* char	outstr[maxch]		# o: column name, units, or print format */
/* int	maxch			# i: maximum length of output string */

void c_tbpset(IRAFPointer tp, int setwhat, int value);
/* pointer tp			# i: pointer to table descriptor */
/* int	setwhat			# i: specifies what parameter is to be set */
/* int	value			# i: the value that is to be assigned */

int c_tbpsta(IRAFPointer tp, int param);
/* pointer tp			# i: pointer to table descriptor */
/* int	param			# i: the parameter to be determined. */

IRAFPointer c_tbcdes(IRAFPointer tp, IRAFPointer cp);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to column descriptor */

void c_tbsirow(IRAFPointer tp, int selrow, int *rownum);
/* pointer tp		# i: pointer to table descriptor */
/* int	selrow		# i: row number (or selected row number) */
/* int	rownum		# o: actual row number */

void c_tbegtb(IRAFPointer tp, IRAFPointer cptr, int rownum, Bool *buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* bool	buffer			# o: buffer for value to be gotten */
/* int	rownum			# actual row number */

void c_tbegtd(IRAFPointer tp, IRAFPointer cptr, int rownum, double *buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* double	buffer			# o: buffer for value to be gotten */
/* int	rownum			# actual row number */

void c_tbegtr(IRAFPointer tp, IRAFPointer cptr, int rownum, float *buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* real	buffer			# o: buffer for value to be gotten */
/* int	rownum			# actual row number */

void c_tbegti(IRAFPointer tp, IRAFPointer cptr, int rownum, int *buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* int	buffer			# o: buffer for value to be gotten */
/* int	rownum			# actual row number */

void c_tbegts(IRAFPointer tp, IRAFPointer cptr, int rownum, short *buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* short	buffer			# o: buffer for value to be gotten */
/* int	rownum			# actual row number */

void c_tbegtt(IRAFPointer tp, IRAFPointer cptr, int rownum, char *buffer, int maxch);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* char	buffer[ARB]		# o: buffer for value to be gotten */
/* int	maxch			# i: max number of char in output string */
/* int	rownum			# actual row number */

void c_tbeptb(IRAFPointer tp, IRAFPointer cptr, int rownum, Bool buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* bool	buffer			# i: value to be put */
/* int	rownum			# actual row number */

void c_tbeptd(IRAFPointer tp, IRAFPointer cptr, int rownum, double buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* double	buffer			# i: value to be put */
/* int	rownum			# actual row number */

void c_tbeptr(IRAFPointer tp, IRAFPointer cptr, int rownum, float buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* real	buffer			# i: value to be put */
/* int	rownum			# actual row number */

void c_tbepti(IRAFPointer tp, IRAFPointer cptr, int rownum, int buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* int	buffer			# i: value to be put */
/* int	rownum			# actual row number */

void c_tbepts(IRAFPointer tp, IRAFPointer cptr, int rownum, short buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* short	buffer			# i: value to be put */
/* int	rownum			# actual row number */

void c_tbeptt(IRAFPointer tp, IRAFPointer cptr, int rownum, char *buffer);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cptr			# i: pointer to column descriptor */
/* char	buffer[ARB]		# i: value to be put */
/* int	rownum			# actual row number */

int c_tbagtd(IRAFPointer tp, IRAFPointer cp, int row, double *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* double	buffer[ARB]	# o: values read from table */
/* int	first		# i: number of first array element to read */
/* int	nelem		# i: maximum number of elements to read */
/* int	row		# actual row number */

int c_tbagtr(IRAFPointer tp, IRAFPointer cp, int row, float *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* real	buffer[ARB]	# o: values read from table */
/* int	first		# i: number of first array element to read */
/* int	nelem		# i: maximum number of elements to read */
/* int	row		# actual row number */

int c_tbagti(IRAFPointer tp, IRAFPointer cp, int row, int *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* int	buffer[ARB]	# o: values read from table */
/* int	first		# i: number of first array element to read */
/* int	nelem		# i: maximum number of elements to read */
/* int	row		# actual row number */

int c_tbagts(IRAFPointer tp, IRAFPointer cp, int row, short *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* short	buffer[ARB]	# o: values read from table */
/* int	first		# i: number of first array element to read */
/* int	nelem		# i: maximum number of elements to read */
/* int	row		# actual row number */

int c_tbagtb(IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* bool	buffer[ARB]	# o: values read from table */
/* int	first		# i: number of first array element to read */
/* int	nelem		# i: maximum number of elements to read */
/* int	row		# actual row number */

int c_tbagtt(IRAFPointer tp, IRAFPointer cp, int row, char **cbuf, int maxch, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* char	cbuf[maxch,ARB]	# o: values read from table */
/* int	maxch		# i: size of first dimension of cbuf */
/* int	first		# i: number of first array element to read */
/* int	nelem		# i: maximum number of elements to read */
/* int	row		# actual row number */

void c_tbaptd(IRAFPointer tp, IRAFPointer cp, int row, double *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* double	buffer[ARB]	# i: values to write to table */
/* int	first		# i: number of first array element to write */
/* int	nelem		# i: number of elements to write */
/* int	row		# actual row number */

void c_tbaptr(IRAFPointer tp, IRAFPointer cp, int row, float *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* real	buffer[ARB]	# i: values to write to table */
/* int	first		# i: number of first array element to write */
/* int	nelem		# i: number of elements to write */
/* int	row		# actual row number */

void c_tbapti(IRAFPointer tp, IRAFPointer cp, int row, int *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* int	buffer[ARB]	# i: values to write to table */
/* int	first		# i: number of first array element to write */
/* int	nelem		# i: number of elements to write */
/* int	row		# actual row number */

void c_tbapts(IRAFPointer tp, IRAFPointer cp, int row, short *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* short	buffer[ARB]	# i: values to write to table */
/* int	first		# i: number of first array element to write */
/* int	nelem		# i: number of elements to write */
/* int	row		# actual row number */

void c_tbaptb(IRAFPointer tp, IRAFPointer cp, int row, Bool *buffer, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* bool	buffer[ARB]	# i: values to write to table */
/* int	first		# i: number of first array element to write */
/* int	nelem		# i: number of elements to write */
/* int	row		# actual row number */

void c_tbaptt(IRAFPointer tp, IRAFPointer cp, int row, char **cbuf, int maxch, int first, int nelem);
/* pointer tp		# i: pointer to table struct */
/* pointer cp		# i: pointer to column struct */
/* char	cbuf[maxch,ARB]	# i: values to write to table */
/* int	maxch		# i: size of first dimension of cbuf */
/* int	first		# i: number of first array element to write */
/* int	nelem		# i: number of elements to write */
/* int	row		# actual row number */

void c_tbrgtb(IRAFPointer tp, IRAFPointer *cp, Bool *buffer, Bool *nullflag, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* bool	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	numcols			# i: number of columns from which to get values */
/* int	rownum			# actual row number */

void c_tbrgtd(IRAFPointer tp, IRAFPointer *cp, double *buffer, Bool *nullflag, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* double	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	numcols			# i: number of columns from which to get values */
/* int	rownum			# actual row number */

void c_tbrgtr(IRAFPointer tp, IRAFPointer *cp, float *buffer, Bool *nullflag, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* real	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	numcols			# i: number of columns from which to get values */
/* int	rownum			# actual row number */

void c_tbrgti(IRAFPointer tp, IRAFPointer *cp, int *buffer, Bool *nullflag, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* int	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	numcols			# i: number of columns from which to get values */
/* int	rownum			# actual row number */

void c_tbrgts(IRAFPointer tp, IRAFPointer *cp, short *buffer, Bool *nullflag, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* short	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	numcols			# i: number of columns from which to get values */
/* int	rownum			# actual row number */

void c_tbrgtt(IRAFPointer tp, IRAFPointer *cp, char **buffer, Bool *nullflag, int lenstr, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* char	buffer[lenstr,ARB]	# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	lenstr			# i: length of each string in array buffer */
/* int	numcols			# i: number of columns from which to get values */
/* int	rownum			# actual row number */

void c_tbrptb(IRAFPointer tp, IRAFPointer *cp, Bool *buffer, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* bool	buffer[ARB]		# i: array of values to be put into table */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

void c_tbrptd(IRAFPointer tp, IRAFPointer *cp, double *buffer, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* double	buffer[ARB]		# i: array of values to be put into table */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

void c_tbrptr(IRAFPointer tp, IRAFPointer *cp, float *buffer, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* real	buffer[ARB]		# i: array of values to be put into table */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

void c_tbrpti(IRAFPointer tp, IRAFPointer *cp, int *buffer, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* int	buffer[ARB]		# i: array of values to be put into table */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

void c_tbrpts(IRAFPointer tp, IRAFPointer *cp, short *buffer, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* short	buffer[ARB]		# i: array of values to be put into table */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

void c_tbrptt(IRAFPointer tp, IRAFPointer *cp, char **buffer, int lenstr, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp[ARB]			# i: array of pointers to column descriptors */
/* char	buffer[lenstr,ARB]	# i: array of values to be put into table */
/* int	lenstr			# i: length of each string in array buffer */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

void c_tbcgtb(IRAFPointer tp, IRAFPointer cp, Bool *buffer, Bool *nullflag, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* bool	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	firstrow, lastrow	# actual range of row numbers */

void c_tbcgtd(IRAFPointer tp, IRAFPointer cp, double *buffer, Bool *nullflag, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* double	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	firstrow, lastrow	# actual range of row numbers */

void c_tbcgtr(IRAFPointer tp, IRAFPointer cp, float *buffer, Bool *nullflag, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* real	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	firstrow, lastrow	# actual range of row numbers */

void c_tbcgti(IRAFPointer tp, IRAFPointer cp, int *buffer, Bool *nullflag, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* int	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	firstrow, lastrow	# actual range of row numbers */

void c_tbcgts(IRAFPointer tp, IRAFPointer cp, short *buffer, Bool *nullflag, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* short	buffer[ARB]		# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	firstrow, lastrow	# actual range of row numbers */

void c_tbcgtt(IRAFPointer tp, IRAFPointer cp, char **buffer, Bool *nullflag, int lenstr, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* char	buffer[lenstr,ARB]	# o: buffer for values */
/* bool	nullflag[ARB]		# o: true if element is undefined in table */
/* int	lenstr			# i: length of each element of buffer */
/* int	firstrow, lastrow	# actual range of row numbers */

void c_tbcptb(IRAFPointer tp, IRAFPointer cp, Bool *buffer, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to descriptor of the column */
/* bool	buffer[ARB]	# i: array of values to be put into column */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbcptd(IRAFPointer tp, IRAFPointer cp, double *buffer, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to descriptor of the column */
/* double	buffer[ARB]	# i: array of values to be put into column */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbcptr(IRAFPointer tp, IRAFPointer cp, float *buffer, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to descriptor of the column */
/* real	buffer[ARB]	# i: array of values to be put into column */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbcpti(IRAFPointer tp, IRAFPointer cp, int *buffer, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to descriptor of the column */
/* int	buffer[ARB]	# i: array of values to be put into column */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbcpts(IRAFPointer tp, IRAFPointer cp, short *buffer, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to descriptor of the column */
/* short	buffer[ARB]	# i: array of values to be put into column */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbcptt(IRAFPointer tp, IRAFPointer cp, char **buffer, int lenstr, int firstrow, int lastrow);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to descriptor of the column */
/* char	buffer[lenstr,ARB]	# i: array of values to be put into column */
/* int	lenstr			# i: number of char in each element of buffer */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbrudf(IRAFPointer tp, IRAFPointer *colptr, int numcols, int rownum);
/* pointer tp			# i: pointer to table descriptor */
/* pointer colptr[ARB]		# i: array of pointers to column descriptors */
/* int	numcols			# i: number of columns */
/* int	rownum			# actual row number */

Bool c_tbhgtb(IRAFPointer tp, char *keyword);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[ARB]		# i: name of parameter to get */

double c_tbhgtd(IRAFPointer tp, char *keyword);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[ARB]		# i: name of parameter to get */

int c_tbhgti(IRAFPointer tp, char *keyword);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[ARB]		# i: name of parameter to get */

float c_tbhgtr(IRAFPointer tp, char *keyword);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[ARB]		# i: name of parameter to get */

void c_tbhgtt(IRAFPointer tp, char *keyword, char *text, int maxch);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[ARB]		# i: name of parameter */
/* char	text[ARB]		# o: value of parameter */
/* int	maxch			# i: maximum number of characters to get */

void c_tbhadb(IRAFPointer tp, char *keyword, Bool value);
/* pointer tp			# i: Pointer to table descriptor */
/* bool	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhadd(IRAFPointer tp, char *keyword, double value);
/* pointer tp			# i: Pointer to table descriptor */
/* double	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhadi(IRAFPointer tp, char *keyword, int value);
/* pointer tp			# i: Pointer to table descriptor */
/* int	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhadr(IRAFPointer tp, char *keyword, float value);
/* pointer tp			# i: Pointer to table descriptor */
/* real	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhadt(IRAFPointer tp, char *keyword, char *text);
/* pointer tp			# i: Pointer to table descriptor */
/* char	keyword[ARB]		# i: Name of parameter */
/* char	text[ARB]		# i: Value of parameter */

void c_tbhptb(IRAFPointer tp, char *keyword, Bool value);
/* pointer tp			# i: Pointer to table descriptor */
/* bool	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhptd(IRAFPointer tp, char *keyword, double value);
/* pointer tp			# i: Pointer to table descriptor */
/* double	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhpti(IRAFPointer tp, char *keyword, int value);
/* pointer tp			# i: Pointer to table descriptor */
/* int	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhptr(IRAFPointer tp, char *keyword, float value);
/* pointer tp			# i: Pointer to table descriptor */
/* real	value			# i: Value of parameter */
/* char	keyword[ARB]		# i: Name of parameter */

void c_tbhptt(IRAFPointer tp, char *keyword, char *text);
/* pointer tp			# i: Pointer to table descriptor */
/* char	keyword[ARB]		# i: Name of parameter */
/* char	text[ARB]		# i: Value of parameter */

void c_tbhcal(IRAFPointer itp, IRAFPointer otp);
/* pointer itp			# i: pointer to descriptor of input table */
/* pointer otp			# i: pointer to descriptor of output table */

void c_tbhgnp(IRAFPointer tp, int parnum, char *keyword, int *dtype, char *str);
/* pointer tp			# i: pointer to table descriptor */
/* int	parnum			# i: number of the parameter to be gotten */
/* char	keyword[SZ_KEYWORD]	# o: keyword for the parameter */
/* int	dtype			# o: data type (TY_CHAR, etc) */
/* char	str[SZ_PARREC]		# o: string containing the value of the param. */

void c_tbhgcm(IRAFPointer tp, char *keyword, char *comment, int maxch);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[SZ_KEYWORD]	# i: keyword to be found */
/* char	comment[ARB]		# o: comment string for keyword */
/* int	maxch			# i: max size of comment */

void c_tbhpcm(IRAFPointer tp, char *keyword, char *comment);
/* pointer tp			# i: pointer to table descriptor */
/* char	keyword[SZ_KEYWORD]	# i: keyword to be found */
/* char	comment[ARB]		# i: comment string for keyword */

void c_tbfpri(char *intable, char *outtable, int *copied);
/* char	intable[ARB]	# i: name of FITS file to be copied */
/* char	outtable[ARB]	# i: name of new FITS file */
/* int	copied		# o: YES if input header was copied to output */

void c_tbtcpy(char *intable, char *outtable);
/* char	intable[ARB]		# i: name of table to be copied to outtable */
/* char	outtable[ARB]		# i: name of new table */

void c_tbtdel(char *table);
/* char	table[ARB]	# i: name of table to be deleted */

void c_tbtren(char *intable, char *outtable);
/* char	intable[ARB]	# i: name of table to be renamed to outtable */
/* char	outtable[ARB]	# i: new name of table */

int c_tbtacc(char *tablename);
/* char	tablename[ARB]		# i: the table name */

void c_tbtext(char *inname, char *outname, int maxch);
/* char	inname[ARB]		# i: table name, possibly without extension */
/* char	outname[ARB]		# o: table name, including extension */
/* int	maxch			# i: max number of char in inname or outname */

void c_tbtnam(IRAFPointer tp, char *tblname, int maxch);
/* pointer tp			# i: pointer to table descriptor */
/* char	tblname[ARB]		# o: the name of the table */
/* int	maxch			# i: maximum number of characters in name */

void c_tbtflu(IRAFPointer tp);
/* pointer tp			# i: pointer to table descriptor */

void c_tbtchs(IRAFPointer tp, int maxpar, int maxcols, int rowlen, int allrows);
/* pointer tp			# i: pointer to table descriptor */
/* int	maxpar			# i: new value for max number of header keywords */
/* int	maxcols			# i: new value for maximum number of columns */
/* int	rowlen			# i: new value for row length */
/* int	allrows			# i: new value of allocated number of rows */

void c_tbrcpy(IRAFPointer itp, IRAFPointer otp, int iselrow, int oselrow);
/* pointer itp			# i: pointer to descriptor of input table */
/* pointer otp			# i: pointer to descriptor of output table */
/* int	iselrow			# i: row number (selected row) in input table */
/* int	oselrow			# i: row number (selected row) in output table */

void c_tbrcsc(IRAFPointer itp, IRAFPointer otp, IRAFPointer *icp, IRAFPointer *ocp, int irow, int orow, int ncols);
/* pointer itp		# i: pointer to descriptor of input table */
/* pointer otp		# i: pointer to descriptor of output table */
/* pointer icp[ncols]	# i: array of pointers for input columns */
/* pointer ocp[ncols]	# i: array of pointers for output columns */
/* int	irow		# i: row number in input table */
/* int	orow		# i: row number in output table */
/* int	ncols		# i: number of columns to be copied */

void c_tbrswp(IRAFPointer tp, int row1, int row2);
/* pointer tp		# i: pointer to table descriptor */
/* int	row1, row2		# actual row numbers corresponding to selrow1,2 */

void c_tbtsrt(IRAFPointer tp, int numcols, IRAFPointer *colptr, Bool fold, int nindex, int *index);
/* pointer	tp		#  i: Table descriptor */
/* int	numcols		#  i: Number of columns to sort on */
/* pointer	colptr[ARB]	#  i: Array of column descriptors */
/* bool	fold		#  i: Fold upper and lower case when sorting */
/* int	nindex		#  i: Number of rows */
/* int	index[ARB]	# io: Array of row indices in sorted order */

void c_tbrdel(IRAFPointer tp, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* int	firstrow	# i: first row to be deleted */
/* int	lastrow		# i: last row to be deleted */

void c_tbrnll(IRAFPointer tp, int firstrow, int lastrow);
/* pointer tp		# i: pointer to table descriptor */
/* int	firstrow, lastrow	# range of actual row numbers */

void c_tbcnam(IRAFPointer tp, IRAFPointer cp, char *colname);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to a column descriptor */
/* char	colname[ARB]		# i: column name */

void c_tbcfmt(IRAFPointer tp, IRAFPointer cp, char *colfmt);
/* pointer tp			# i: pointer to table descriptor */
/* pointer cp			# i: pointer to a column descriptor */
/* char	colfmt[ARB]		# i: print format for column */

void c_tbcnit(IRAFPointer tp, IRAFPointer cp, char *colunits);
/* pointer tp		# i: pointer to table descriptor */
/* pointer cp		# i: pointer to a column descriptor */
/* char	colunits[ARB]	# i: new value of units for column */

IRAFPointer c_tbcnum(IRAFPointer tp, int colnum);
/* pointer tp			# i: pointer to table descriptor */
/* int	colnum			# i: column number (not pointer) */

int c_tbparse(char *tablename, char *fname, char *extname, int maxch, int *hdu);
/* char	tablename[ARB]	# i: name as specified by user */
/* char	fname[ARB]	# o: name of file containing table */
/* char	extname[ARB]	# o: CDF name, or null if none */
/* int	maxch		# i: size of fname and extname strings */
/* int	hdu		# o: HDU number for FITS file, or -1 if none */

