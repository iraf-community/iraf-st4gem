# STR_COMPARE -- Comparison function for two integers, possibly INDEF
#
# B.Simon	29-Jun-88	First Code

int procedure str_compare (obj, size, i, j)

pointer	obj		# i: Pointer to a string array
int	size		# i: Length of the string
int	i		# i: Index of first element to be compared
int	j		# i: Index of second element to be compared
#--
int	order

bool	streq(), strgt()

begin
	# Indefinite values are defined to be greater than
	# those that are not indefinite

	if (streq (Memc[obj+(i-1)*(size+1)], Memc[obj+(j-1)*(size+1)]))
	    order = 0
	else if (Memc[obj+(i-1)*(size+1)] == EOS)
	    order = 1
	else if (Memc[obj+(j-1)*(size+1)] == EOS)
	    order = -1
	else if (strgt (Memc[obj+(i-1)*(size+1)], Memc[obj+(j-1)*(size+1)]))
	    order = 1
	else 
	    order = -1

	return (order)
end
