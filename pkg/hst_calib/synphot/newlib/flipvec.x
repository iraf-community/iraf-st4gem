# FLIPVEC -- Reverse the order of the elements in a vector

procedure flipvec (n, x)

int	n		# i: vector length
real	x[ARB]		# u: vector
#--
int	i, j
real	temp

begin
	i = 1
	j = n

	while (i < j) {
	    temp = x[i]
	    x[i] = x[j]
	    x[j] = temp

	    i = i + 1
	    j = j - 1
	}

end
