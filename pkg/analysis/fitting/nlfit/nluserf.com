# NLUSERF.COM -- Global variables used by user function.

pointer	lx, lcoef
int	lndata, lncoef
bool	optimize

common	/nlufcom/  lx, lcoef, lndata, lncoef, optimize

