# MERGESORT -- General merge sort for arbitrary objects.  X is an integer array
# indexing the array to be sorted.  The user supplied COMPARE function is used
# to compare objects indexed by X:
# 
# 	-1,0,1 = compare (obj, x1, x2)
# 
# where the value returned by COMPARE has the following significance:
# 
# 	-1	obj[x1]  < obj[x2]
# 	 0	obj[x1] == obj[x2]
# 	 1	obj[x1]  > obj[x2]
# 
# MSORT reorders the elements of the X array, which must be of type integer.
#
# B.Simon	28-Sept-87	First Code
# B.Simon	29-Jun-88	Modified for graph tasks

procedure mergesort (x, obj, size, mx, nx, compare)

int	x[ARB]		# io: Index array to be sorted
pointer	obj		#  i: Pointer to array used in comparisons
int	size		#  i: Length of the array element
int	mx		#  i: Length of array x (Must be >= 2 * nx)
int	nx		#  i: Number of elements to be sorted
extern	compare()	#  i: Function to be called to compare elements
#--
bool	up
int	ielem, jelem, kelem, melem
int	runlen, ilen, jlen

string	toosmall  "mergesort: index array too small"
int	compare()

begin
	if (2 * nx > mx)
	    call printerr_int (toosmall, mx)

	# Merging two sorted runs creates a new sorted run twice the length
	# of the original run. Continue this process until the sorted run
	# length is equal to the array length.

	up = false
	for (runlen = 1; runlen < nx; runlen = 2 * runlen) {

	    # The runs are stored in one of two halves of the x array.
	    # Set the array pointers according to the half the runs are
	    # located in now.

	    if (! up) {
		ielem = 1
		jelem = runlen + 1
		kelem = mx - nx + 1
		melem = nx
	    } else {
		ielem = mx - nx + 1
		jelem = runlen + ielem
		kelem = 1
		melem = mx
	    }

	    # Loop over each pair of runs in the array

	    while (ielem <= melem) {
		ilen = min (runlen, melem-ielem+1)
		jlen = min (runlen, melem-jelem+1)

		# Merge the pair of runs into the other half of the x array

		while (ilen > 0 && jlen > 0) {
		    if (compare (obj, size, x[ielem], x[jelem]) <= 0) {
			x[kelem] = x[ielem]
			ielem = ielem + 1
			kelem = kelem + 1
			ilen = ilen - 1
		    } else {
			x[kelem] = x[jelem]
			jelem = jelem + 1
			kelem = kelem + 1
			jlen = jlen - 1
		    }
		}

		# Copy the remaining elements from i when j is exhausted

		while (ilen > 0) {
		    x[kelem] = x[ielem]
		    ielem = ielem + 1
		    kelem = kelem + 1
		    ilen = ilen - 1
		}

		# Copy the remaining elements from j when i is exhausted

		while (jlen > 0) {
		    x[kelem] = x[jelem]
		    jelem = jelem + 1
		    kelem = kelem + 1
		    jlen = jlen - 1
		}

		# Set array pointers to next set of runs

		ielem = ielem + runlen
		jelem = jelem + runlen
	    }
	    up = ! up
	}

	# If result is in the upper end of x array, move it to the lower
	# end

	if (up)
	    call amovi (x[mx-nx+1], x[1], nx)

	end
