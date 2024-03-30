# VZERO.COM -- Common block used to track vzero values

real	vzero[3,MAXLIST]	# first, last, step for a list of ranges
real	vzcur			# current value of vzero
int	vzcount			# count of remaining values in current range
int	vzindex			# index to current range
int	vzlast			# last range in list

common	/ vzcom / vzero, vzcur, vzcount, vzindex, vzlast
