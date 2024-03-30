pointer	buf		# Rectangular buffer to be interpolated on
int	nx		# First dimension of buffer

common  / earcom /  buf, nx
