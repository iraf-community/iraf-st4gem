pointer	buf		# Rectangular buffer to be interpolated on
int	nx		# First dimension of buffer
int	ny		# Second dimension of buffer

common  / lincom /  buf, nx, ny
