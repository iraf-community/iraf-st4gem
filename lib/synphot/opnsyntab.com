# OPNSYTAB.COM -- Common block for persistent variables used by opnsyntab

pointer	tabptr[MAXCACHE]    # array of table pointers
long	lastuse[MAXCACHE]   # value of counter last time table was opened
long	nextuse		    # counter for least recently used algorithm
int	initial		    # YES if common block was initialized

common	/syntab/	tabptr, lastuse, nextuse, initial

