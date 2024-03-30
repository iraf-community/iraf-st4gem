# INT_COMPARE -- Comparison function for two integers, possibly INDEF
#
# B.Simon	29-Jun-88	First Code

int procedure int_compare (obj, size, i, j)

pointer	obj		# i: Pointer to an integer array
int	size		# i: (Not used)
int	i		# i: Index of first element to be compared
int	j		# i: Index of second element to be compared
#--
int	order

begin
	# Indefinite values are defined to be greater than
	# those that are not indefinite

	if (Memi[obj+i-1] == Memi[obj+j-1])
	    order = 0
	else if (Memi[obj+i-1] == INDEFI)
	    order = 1
	else if (Memi[obj+j-1] == INDEFI)
	    order = -1
	else if (Memi[obj+i-1] > Memi[obj+j-1])
	    order = 1
	else 
	    order = -1

	return (order)
end
