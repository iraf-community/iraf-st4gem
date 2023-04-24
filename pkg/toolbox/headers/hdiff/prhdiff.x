include <imio.h>
include "hdiff.h"

# PR_HDIFF -- Print differences between two header files
#
# This procedure compares two sorted lists of header keywords. A line is
# printed if a keyword is found in one list and not in the other, or the
# value of the two keywords differs. Since the two lists are sorted, a 
# standard merge algorithm can be used to find the differences.
#
# B.Simon	25-Apr-90	Original

procedure pr_hdiff (im1, im2, buffer1, index1, ndx1, buffer2, index2, ndx2)

pointer	im1		# i: First image descriptor
pointer	im2		# i: Second image descriptor
char	buffer1[ARB]	# i: First buffer of keyword records
int	index1[ARB]	# i: First array of keyword record indices
int	ndx1		# i: Number of records in first buffer
char	buffer2[ARB]	# i: Second buffer of keyword records
int	index2[ARB]	# i: Second array of keyword record indices
int	ndx2		# i: Number of records in second buffer
#--
bool	title
int	idx1, idx2, jdx1, jdx2, order

string	nullstring  ""

int	strncmp(), cmp_record()

begin
	title = true

	# Since it is the index array that is sorted, and not the string
	# buffer, two indices are needed to access the next record

	idx1 = 1
	idx2 = 1
	jdx1 = index1[idx1]
	jdx2 = index2[idx2]

	# Main phase of merge, when neither list of records is empty

	while (idx1 <= ndx1 && idx2 <= ndx2) {
	    order = strncmp (buffer1[jdx1], buffer2[jdx2], SZ_KEYWORD)
	    if (order < 0) {
		if (title) {
		    title = false
		    call pr_title (im1, im2)
		}
		call pr_record (buffer1[jdx1], nullstring)
		idx1 = idx1 + 1
		jdx1 = index1[idx1]
	    } else if (order > 0) {
		if (title) {
		    title = false
		    call pr_title (im1, im2)
		}
		call pr_record (nullstring, buffer2[jdx2])
		idx2 = idx2 + 1
		jdx2 = index2[idx2]
	    } else {
		if (cmp_record (buffer1[jdx1], buffer2[jdx2]) != 0) {
		    if (title) {
			title = false
			call pr_title (im1, im2)
			}
		    call pr_record (buffer1[jdx1], buffer2[jdx2])
		}
		idx1 = idx1 + 1
		jdx1 = index1[idx1]
		idx2 = idx2 + 1
		jdx2 = index2[idx2]		
	    }
	}

	# Print the tail of the longer list. At most, only one of
	# these while loops will be called

	while (idx1 <= ndx1) {
	    if (title) {
		title = false
		call pr_title (im1, im2)
	    }
	    call pr_record (buffer1[jdx1], nullstring)
	    idx1 = idx1 + 1
	    jdx1 = index1[idx1]
	}

	while (idx2 <= ndx2) {
	    if (title) {
		title = false
		call pr_title (im1, im2)
	    }
	    call pr_record (buffer2[jdx2], nullstring)
	    idx2 = idx2 + 1
	    jdx2 = index2[idx1]
	}

end

# PR_TITLE -- Print title containing image names

procedure pr_title (im1, im2)

pointer	im1		# i: First image descriptor
pointer im2		# i: Second image descriptor
#--
int	offset1, offset2
pointer	sp, file1, file2, temp, root1, root2

int	fnldir()

begin
	call smark (sp)
	call salloc (file1, SZ_FNAME, TY_CHAR)
	call salloc (file2, SZ_FNAME, TY_CHAR)
	call salloc (temp, SZ_FNAME, TY_CHAR)

	# Remove the directory name from the root, if present

	call imgcluster (IM_NAME(im1), Memc[file1], SZ_FNAME)
	offset1 = fnldir (Memc[file1], Memc[temp], SZ_FNAME)
	call strcpy (IM_NAME(im1), Memc[file1], SZ_FNAME)
	root1 = file1 + offset1

	call imgcluster (IM_NAME(im2), Memc[file2], SZ_FNAME)
	offset2 = fnldir (Memc[file2], Memc[temp], SZ_FNAME)
	call strcpy (IM_NAME(im2), Memc[file2], SZ_FNAME)
	root2 = file2 + offset2

	call printf ("#\n#%11t%-30s%-30s\n#\n")
	call pargstr (Memc[root1])
	call pargstr (Memc[root2])

	call sfree (sp)
end

# PR_RECORD -- Print the two records that differ

procedure pr_record (record1, record2)

char	record1[ARB]	# i: Keyword record from first image
char	record2[ARB]	# i: Keyword record from second image
#--
int	ic, nc
pointer	sp, keyword, value1, value2

int	ctowrd()

begin
	call smark (sp)
	call salloc (keyword, SZ_KEYWORD, TY_CHAR)
	call salloc (value1, SZ_VALUE, TY_CHAR)
	call salloc (value2, SZ_VALUE, TY_CHAR)

	# One or the other of the records may be a null string, 
	# indicating that this record was not present in the image.

	if (record1[1] != EOS) {
	    call strcpy (record1, Memc[keyword], SZ_KEYWORD)
	} else if (record2[1] != EOS) {
	    call strcpy (record2, Memc[keyword], SZ_KEYWORD)
	} else {
	    call error (1, "pr_record: keyword not found")
	}

	# Get the value field from the record and trim trailing blanks

	if (record1[1] == EOS) {
	    Memc[value1] = EOS
	} else {
	    ic = START_VALUE
	    nc = ctowrd (record1, ic, Memc[value1], SZ_VALUE)
	    for (ic = nc - 1; nc >= 0; ic = ic - 1) {
		if (Memc[value1+ic] != ' ')
		    break
	    }
	    Memc[value1+ic+1] = EOS
	}

	if (record2[1] == EOS) {
	    Memc[value2] = EOS
	} else {
	    ic = START_VALUE
	    nc = ctowrd (record2, ic, Memc[value2], SZ_VALUE)
	    for (ic = nc - 1; nc >= 0; ic = ic - 1) {
		if (Memc[value2+ic] != ' ')
		    break
	    }
	    Memc[value2+ic+1] = EOS
	}

	# Write a line containing the keyword names and the two values,
	# either of which may be a null string

	call printf ("%-10s%-30s%-30s\n")
	call pargstr (Memc[keyword])
	call pargstr (Memc[value1])
	call pargstr (Memc[value2])

	call sfree (sp)
end
