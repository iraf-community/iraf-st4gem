# include <c_iraf.h>

# define CL_PRTYPE	1
# define CL_PCACHE	2

# define PR_CONNECTED 1
# define PR_DETACHED  2
# define PR_HOST      3

# define PR_NOEXIT	0
# define PR_EXIT	1	

void c_clseti(int parameter, int value);
/* int	parameter		# CLIO parameter being queried */
/* int	value			# value of parameter (output) */

int c_clstati(int parameter);
/* int	parameter */

Bool c_clgetb(char *param);
/* char	param[ARB] */

char c_clgetc(char *param);
/* char	param[ARB] */

double c_clgetd(char *param);
/* char	param[ARB] */

int c_clgeti(char *param);
/* char	param[ARB] */

int c_clgetl(char *param);
/* char	param[ARB] */

float c_clgetf(char *param);
/* char	param[ARB] */

float c_clgetr(char *param);
/* char	param[ARB] */

short c_clgets(char *param);
/* char	param[ARB] */

Complex c_clgetx(char *param);
/* char	param[ARB] */

void c_clputb(char *param, Bool bval);
/* char	param[ARB] */
/* bool	bval */

void c_clputc(char *param, char cval);
/* char	param[ARB] */
/* char	cval */

void c_clputd(char *param, double dval);
/* char	param[ARB] */
/* double	dval */

void c_clputi(char *param, int value);
/* char	param[ARB]		# parameter to be set */
/* int	value			# new value */

void c_clputs(char *param, short value);
/* char	param[ARB]		# parameter to be set */
/* short	value			# new value */

void c_clputl(char *param, int lval);
/* char	param[ARB] */
/* long	lval */

void c_clputf(char *param, float rval);
/* char	param[ARB] */
/* real	rval */

void c_clputr(char *param, float rval);
/* char	param[ARB] */
/* real	rval */

void c_clputx(char *param, Complex xval);
/* char	param[ARB] */
/* complex	xval */

void c_clgstr(char *param, char *outstr, int maxch);
/* char	param[ARB], outstr[maxch] */
/* int	maxch */

void c_clpstr(char *param, char *value);
/* char	param[ARB]		# param name */
/* char	value[ARB]		# new value */

int c_clglpb(char *param, Bool *bval);
/* char	param[ARB] */
/* bool	bval */

int c_clglpc(char *param, char *cval);
/* char	param[ARB] */
/* char	cval */

int c_clglpd(char *param, double *dval);
/* char	param[ARB] */
/* double	dval */

int c_clglpi(char *param, int *ival);
/* char	param[ARB] */
/* int	ival, stat, clglpd() */

int c_clglpl(char *param, int *lval);
/* char	param[ARB] */
/* long	lval */

int c_clglpf(char *param, float *rval);
/* char	param[ARB] */
/* real	rval */

int c_clglpr(char *param, float *rval);
/* char	param[ARB] */
/* real	rval */

int c_clglps(char *param, short *sval);
/* char	param[ARB] */
/* short	sval */

int c_clglpx(char *param, Complex *xval);
/* char	param[ARB] */
/* complex	xval */

int c_clglstr(char *param, char *outstr, int maxch);
/* char	param[ARB], outstr[maxch] */
/* int	maxch */

int c_clgcur(char *param, float *wx, float *wy, int *wcs, int *key, char *strval, int maxch);
/* char	param[ARB]		# parameter to be read */
/* real	wx, wy			# cursor coordinates */
/* int	wcs			# wcs to which coordinates belong */
/* int	key			# keystroke value of cursor event */
/* char	strval[ARB]		# string value, if any */
/* int	maxch */

int c_clgkey(char *param, int *key, char *strval, int maxch);
/* char	param[ARB]		# parameter to be read */
/* int	key			# keystroke value of cursor event */
/* char	strval[ARB]		# string value, if any */
/* int	maxch */

int c_clgwrd(char *param, char *keyword, int maxchar, char *dictionary);
/* char	param[ARB]		# CL parameter string */
/* char	keyword[ARB]		# String matched in dictionary */
/* int	maxchar			# Maximum size of str */
/* char	dictionary[ARB]		# Dictionary string */

IRAFPointer c_clopset(char *pset);
/* char	pset[ARB]		# pset name (name of CL pset parameter) */

void c_clcpset(IRAFPointer pp);
/* pointer	pp		# pset descriptor */

Bool c_clgpsetb(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

char c_clgpsetc(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

double c_clgpsetd(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

int c_clgpseti(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

int c_clgpsetl(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

float c_clgpsetf(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

float c_clgpsetr(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

short c_clgpsets(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

Complex c_clgpsetx(IRAFPointer pp, char *parname);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */

void c_clppsetb(IRAFPointer pp, char *parname, Bool bval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* bool	bval			# new value of parameter */

void c_clppsetc(IRAFPointer pp, char *parname, char cval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* char	cval			# new value of parameter */

void c_clppsetd(IRAFPointer pp, char *parname, double dval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* double	dval			# new value of parameter */

void c_clppseti(IRAFPointer pp, char *parname, int ival);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* int	ival			# new value of parameter */

void c_clppsetl(IRAFPointer pp, char *parname, int lval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* long	lval			# new value of parameter */

void c_clppsetf(IRAFPointer pp, char *parname, float rval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* real	rval			# new value of parameter */

void c_clppsetr(IRAFPointer pp, char *parname, float rval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* real	rval			# new value of parameter */

void c_clppsets(IRAFPointer pp, char *parname, short sval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* short	sval			# new value of parameter */

void c_clppsetx(IRAFPointer pp, char *parname, Complex xval);
/* pointer	pp			# pset descriptor */
/* char	parname[ARB]		# parameter name */
/* complex	xval			# new value of parameter */

void c_clgpseta(IRAFPointer pp, char *pname, char *outstr, int maxch);
/* pointer	pp			# pset descriptor */
/* char	pname[ARB]		# parameter name */
/* char	outstr[maxch]		# output string */
/* int	maxch			# max chars out */

void c_clppseta(IRAFPointer pp, char *pname, char *sval);
/* pointer	pp			# pset descriptor */
/* char	pname[ARB]		# parameter name */
/* char	sval[ARB]		# string value of parameter */

void c_cllpset(IRAFPointer pp, int fd, char *format);
/* pointer	pp			#I pset descriptor */
/* int	fd			#I output file */
/* char	format[ARB]		#I format, one %s each for param, value */

void c_clepset(IRAFPointer pp);
/* pointer	pp			#I pset descriptor */

void c_clcmd(char *cmd);
/* char	cmd[ARB] */

void c_clcmdw(char *cmd);
/* char	cmd[ARB] */

