# IOTA -- Initialize an array to 1 ... n

procedure iota (a, n)

int	a[ARB]		# o: array to be initiialized
int	n		# i: length of array
#--
int	i

begin
	do i = 1, n
	    a[i] = i
end
