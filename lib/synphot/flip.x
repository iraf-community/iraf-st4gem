# FLIP -- Reverse the order of elements in an array

procedure flip (n, a)

int	n		# i: Number of elements in array
real	a[ARB]		# u: Array to be flipped
#--
int	i, j
real	temp

begin
	i = 1
	j = n

	while (i < j) {
	    temp = a[i]
	    a[i] = a[j]
	    a[j] = temp

	    i = i + 1
	    j = j - 1
	}
end
