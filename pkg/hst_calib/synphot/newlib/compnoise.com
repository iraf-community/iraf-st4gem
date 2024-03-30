# COMPNOISE.COM -- Global variables used by noise expression compiler

pointer	icode		# Pointer to next available code
int	ncode		# Number of available codes left
int	readout		# Number of readouts
real	time		# Exposure time

common	/noyzcom/  icode, ncode, readout, time
