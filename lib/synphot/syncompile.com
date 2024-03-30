# SYNCOMPILE.COM -- Global variables used by synphot expression compiler

bool	strict		# Use strict rule for determining end of token
pointer	icode		# Pointer to next available code
int	ncode		# Number of available codes left

common	/syncom/  strict, icode, ncode
