#include "libsynphot.h"     #convention precludes including thi
define  MAXPARAM    3
define  SZ_KEYWRD       32              # Length of keyword
define  PCH             '#'             # Leading character of parameter

bool procedure getkeyval( key, val, keylist, vallist, nvals)

char	key[SZ_KEYWRD]		# i: The key you're searching for
real	val                     # o: The value returned. -1 => key does not exist
char	keylist[SZ_KEYWRD,ARB]  # i: The list of keys in the "dictionary"
real	vallist[MAXPARAM,ARB]        # i: The list of values in the "dictionary"
int	nvals                   # i: The number of entries in the "dictionary"
                                # return value: true if key found, false otherwise

# local variables
int	i
char    pchar
data	pchar	/ PCH /

# functions
bool	streq()

# First check for null dictionary: ahivi(tnparam)>0. do that outside though.
# also outside, make sure keys are all lower case.

#NOTE discrepancy in assumptions here:
# vallist can theoretically have MAXPARAM values for each key.
# however to pick off the correct one, we'd have to pass in tnparam, which
# we don't want to do. 
# Therefore, we take advantage of the fact that for thermback, there will
# never be more than one parameter (value) per keyword.
#
begin

# See whether the key is in the dictionary
	call strlwr(key)  
	call strcat(pchar,key,SZ_KEYWRD)
	do i = 1, nvals {
#	   call eprintf("dictionary id = ...%s...\n")
#	   call pargstr(keylist[1,i])
	   if ( streq(key, keylist[1,i])) {
		val = vallist[1,i]
		return true
	   }
	}

	return false

end
