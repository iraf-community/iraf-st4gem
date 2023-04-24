include	<imhdr.h>
include "hdiff.h"

# GETKEYBUF -- Create a sorted buffer of keyword values

procedure getkeybuf (im, igroup, negate, keynames, buffer, index, nindex)

pointer	im		# i: Image descriptor
int	igroup		# i: Group number
bool	negate		# i: Include keywords NOT in the keyword list
char	keynames[ARB]	# i: List of keyword names
char	buffer[ARB]	# o: String buffer to hold keyword records
int	index[ARB]	# o: Array of pointers to keyword records
int	nindex		# o: Number of keyword records
#--
bool	copy
int	ibuf, keylen
pointer sp, key, iua

string	excluded  "HISTORY,COMMENT"

int	word_match(), gf_gfind()

begin
	# Allocate dynamic memory to hold the keyword

	call smark (sp)
	call salloc (key, SZ_KEYWORD, TY_CHAR)

	ibuf = 1
	nindex = 0

	iua = IM_USERAREA(im)
	while (Memc[iua] != EOS) {

	    # Read keyword name and trim trailing blanks

	    call strcpy (Memc[iua], Memc[key], SZ_KEYWORD)
	    for (keylen = SZ_KEYWORD; keylen > 0; keylen = keylen - 1) {
		if (Memc[key+keylen-1] != ' ')
		    break
	    }
	    Memc[key+keylen] = EOS

	    # Check to see whether this keyword should be copied to the buffer

	    if (keylen == 0) {
		copy = false
	    } else if (negate) {
		copy = word_match (Memc[key], keynames) == 0 &&
		       word_match (Memc[key], excluded) == 0
	    } else {
		copy = word_match (Memc[key], keynames) > 0 &&
		       word_match (Memc[key], excluded) == 0
	    }

	    if (copy && igroup > 1)
		copy = gf_gfind (im, Memc[key]) > 0

	    # Copy record to output buffer

	    if (copy) {
		nindex = nindex + 1
		index[nindex] = ibuf
		while (Memc[iua] != '\n') {
		    buffer[ibuf] = Memc[iua]
		    ibuf = ibuf + 1
		    iua = iua + 1
		}
		buffer[ibuf] = EOS
		ibuf = ibuf + 1
		iua = iua + 1

	    # Otherwise, skip to next record

	    } else {

		while (Memc[iua] != '\n')
		    iua = iua + 1
		iua = iua + 1
	    }
	}

	# Sort the buffer array

	if (nindex > 1)
	    call strsrt (index, buffer, nindex)
	call sfree (sp)
end
