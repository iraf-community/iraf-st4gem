# include <c_iraf.h>

int c_decode_ranges(char *range_string, int *ranges, int max_ranges, int *nvalues);
/* char	range_string[ARB]	# Range string to be decoded */
/* int	ranges[3, max_ranges]	# Range array */
/* int	max_ranges		# Maximum number of ranges */
/* int	nvalues			# The number of values in the ranges */

int c_get_next_number(int *ranges, int *number);
/* int	ranges[ARB]		# Range array */
/* int	number			# Both input and output parameter */

int c_get_previous_number(int *ranges, int *number);
/* int	ranges[ARB]		# Range array */
/* int	number			# Both input and output parameter */

Bool c_is_in_range(int *ranges, int number);
/* int	ranges[ARB]		# Range array */
/* int	number			# Number to be tested against ranges */

int c_isdirectory(char *vfn, char *pathname, int maxch);
/* char	vfn[ARB]		# name to be tested */
/* char	pathname[ARB]		# receives path of directory */
/* int	maxch			# max chars out */

