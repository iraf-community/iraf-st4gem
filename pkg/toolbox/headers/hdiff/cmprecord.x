include "hdiff.h"

# CMP_RECORD -- Compare two record values and determine their ordering
#
# This procedure determines the ordering of the value fields of two keyword
# records. If the value in the first record is less than the value in the
# second, a negative number is returned. If the value in the first record
# is greater than the value in the second, a positive number is returned.
# If they are equal, a zero is returned. Trailing blanks in quoted fields
# are not considered significant in the comparison.
#
# B.Simon	25-Apr-90	Original

int procedure cmp_record (record1, record2)

char	record1[ARB]	# i: First record
char	record2[ARB]	# i: Second record
#--
char	delim1, delim2
int	ic, jc, order

begin
	# Skip over leading blanks and get the value delimeter
	# (either a blank or a single quote)

	for (ic = START_VALUE; record1[ic] == ' '; ic = ic + 1)
	    ;

	if (record1[ic] == '\'') {
	    ic = ic + 1
	    delim1 = '\''
	} else {
	    delim1 = ' '
	}

	for (jc = START_VALUE; record2[jc] == ' '; jc = jc + 1)
	    ;

	if (record2[jc] == '\'') {
	    jc = jc + 1
	    delim2 = '\''
	} else {
	    delim2 = ' '
	}

	# Find the first characters which do not compare
	# or fall off the end of the value

	while (record1[ic] == record2[jc]) {
	    if (record1[ic] == EOS || record1[ic] == delim1)
		break
	    ic = ic + 1
	    jc = jc + 1
	}

	# The two value strings are equal (order == 0) if:
	# Both strings have reached the closing delimeter or
	# One string has reached the closing delimeter and the
	# remaining characters in the other string are blank.
	# Blanks can only be part of a value if the delimeter
	# is a single quote.

	if (record1[ic] == delim1) {
	    if (delim2 == '\'') {
		while (record2[jc] == ' ')
		    jc = jc + 1
	    }
	    if (record2[jc] == delim2) {
		order = 0
	    } else {
		order = - record2[jc]
	    }			
	} else if (record2[jc] == delim2) {
	    if (delim1 == '\'') {
		while (record1[ic] == ' ')
		    ic = ic + 1
	    }
	    if (record1[ic] == delim1) {
		order = 0
	    } else {
		order = record1[ic]
	    }    	
	} else {
	   order = record1[ic] - record2[jc]
	}

	return (order)
end
