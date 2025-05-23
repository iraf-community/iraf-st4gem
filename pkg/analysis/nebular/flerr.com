# flerr.com	--	Include file for passing error message parameters	

# Generic common block for passing error message parameters in 
# FIVEL numerical routines:

real	val1, val2, val3, val4		# error parameters

common	/ferr/	val1, val2, val3, val4

