define	SZ_HSHSTRUCT	4

define	HSH_LEN		Memi[$1]	# hash table length
define	HSH_NXTSTR	Memi[$1+1]	# next string location in buffer
define	HSH_INDEX	Memi[$1+2]	# index into name buffer
define	HSH_BUFFER	Memi[$1+3]	# pointer to name buffer

define	HSH_PTR		Memi[HSH_INDEX($1)+($2)]
define	HSH_STR		Memc[HSH_BUFFER($1)+HSH_PTR($1,$2)]

# These procedures save the names of observation mode keywords to make sure
# that no keyword is printed more than once.
#
# B.Simon	10-Jan-94	original

# ADD_HASH -- Add name to hash table

bool procedure add_hash (hash, name)

pointer	hash		# i: Hash table descriptor
char	name[ARB]	# i: Name to find in table
#--
bool	found
int	ic, key, index

bool	streq()
int	gstrcpy()

begin
	# Calculate hash key

	key = 0
	for (ic =1; name[ic] != EOS; ic = ic + 1)
	    key = key + name[ic]

	# Search thru table until the name or an empty slot is found

	found = false
	index = mod (key, HSH_LEN(hash))
	while (HSH_PTR(hash,index) != NULL) {
	    if (streq (HSH_STR(hash,index), name)) {
		found = true
		break
	    }

	    index = mod (index+1, HSH_LEN(hash))
	}

	# If the name was not found, add it to the table

	if (! found) {
	    HSH_PTR(hash,index) = HSH_NXTSTR(hash)
	    HSH_NXTSTR(hash) = gstrcpy (name, HSH_STR(hash,index), ARB) +
			       HSH_NXTSTR(hash) + 1
	}

	# Return true if the name was added
	# or false if the name was already in the table

	return (! found)
end

# FREE_HASH -- Free memory used by hash table

procedure free_hash (hash)

pointer	hash		# i: Hash table descriptor
#--

begin
	call mfree (HSH_INDEX(hash), TY_INT)
	call mfree (HSH_BUFFER(hash), TY_CHAR)
	call mfree (hash, TY_INT)
end

# OPN_HASH - Allocate new hash table

pointer procedure opn_hash (maxname, maxch)

int	maxname		# i: Maximum number of names in hash table
int	maxch		# i: Maximum length of each name
#--
int	indexlen, buflen
pointer	hash

begin
	# Calculate size of arrays

	buflen = maxname * (maxch + 1) + 1
	indexlen = 2 * maxname

	# Allocate memory

	call malloc (hash, SZ_HSHSTRUCT, TY_INT)
	call malloc (HSH_BUFFER(hash), buflen, TY_CHAR)
	call malloc (HSH_INDEX(hash), indexlen, TY_INT)

	# Initialize hash table

	HSH_NXTSTR(hash) = 1
	HSH_LEN(hash) = indexlen

	call amovki (NULL, HSH_PTR(hash,0), indexlen)
	HSH_STR(hash,0) = EOS

	return (hash)
end


