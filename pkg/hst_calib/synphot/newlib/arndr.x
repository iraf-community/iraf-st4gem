#* HISTORY *
#* B.Simon	30-Mar-95	original

# ARNDR -- Round a real number to the nearest whole number

procedure arndr (a, b, n)

real	a[ARB]		# i: input vector
real	b[ARB]		# o: rounded output vector
int	n		# i: length of vectors
#--
int	i

begin
	do i = 1, n
	    b[i] = aint(a[i])
end
