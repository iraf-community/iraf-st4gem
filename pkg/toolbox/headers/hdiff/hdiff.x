include <imio.h>
include <error.h>
include	"hdiff.h"

# HDIFF -- Show the differences between two header files
#
# B.Simon	23-Apr-90	Original

procedure hdiff ()

#--
pointer	input1		# First image template
pointer	input2		# Second image template
pointer	keywords	# List of keywords to compare
bool	nofile		# Ignore error if files do not exist

bool	negate, bjunk
int	nlist, fd, ic, nc, ilist, junk
int	ngroup1, ngroup2, igroup, len1, len2, ndx1, ndx2
pointer	sp, cluster1, cluster2, image1, image2, keynames
pointer	im1, im2, tp1, tp2, buffer1, buffer2, index1, index2

string	badimcount "The number of images in both lists are not the same"
string	badgcount  "The number of groups in both iamges are not the same"

bool	clgetb(), tp_fetch()
int	imtlen(), open(), getline(), imtgetim(), word_count()
pointer	imtopenp(), tp_open(), immap()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (cluster1, SZ_FNAME, TY_CHAR)
	call salloc (cluster2, SZ_FNAME, TY_CHAR)
	call salloc (image1, SZ_FNAME, TY_CHAR)
	call salloc (image2, SZ_FNAME, TY_CHAR)
	call salloc (keywords, SZ_FNAME, TY_CHAR)
	call salloc (keynames, SZ_COMMAND, TY_CHAR)

	# Read the task parameters

	input1 = imtopenp ("input1")
	input2 = imtopenp ("input2")
	call clgstr ("keywords", Memc[keywords], SZ_FNAME)
	nofile = clgetb ("nofile")

	# Make sure that the number of images is the same in both lists

	nlist = imtlen (input1)
	if (imtlen (input2) != nlist)
	    call error (1, badimcount)

	# Create a list of keyword names

 	if (Memc[keywords] != '@') {
	    call strcpy (Memc[keywords], Memc[keynames], SZ_COMMAND)
	} else {
	    fd = open (Memc[keywords+1], READ_ONLY, TEXT_FILE)
	    for (ic = 0; ic < SZ_COMMAND; ic = ic + nc) {
		nc = getline (fd, Memc[keynames+ic])
		if (nc <= 0)
		    break
	    }
	    Memc[keynames+ic] = EOS
	    call close (fd)
	}

	call strupr (Memc[keynames])
	if (Memc[keynames] == '~') {
	    negate = true
	    keynames = keynames + 1
	} else if (word_count (Memc[keynames]) == 0) {
	    negate = true
	    Memc[keynames] = EOS
	} else {
	    negate = false
	}

	# Compare each image

	do ilist = 1, nlist {
	    junk = imtgetim (input1, Memc[cluster1], SZ_FNAME)
	    junk = imtgetim (input2, Memc[cluster2], SZ_FNAME)

	    # Compare each group in each image

	    iferr {
		tp1 = tp_open (Memc[cluster1], 0, ngroup1)
		tp2 = tp_open (Memc[cluster2], 0, ngroup2)
	    } then if (nofile) {
		next
	    } else {
		call erract (EA_ERROR)
	    }

	    if (ngroup1 != ngroup2)
		call error (1, badgcount)
		    
	    do igroup = 1, ngroup1 {
		bjunk = tp_fetch (tp1, Memc[image1])
		bjunk = tp_fetch (tp2, Memc[image2])

		im1 = immap (Memc[image1], READ_ONLY, 0)
		im2 = immap (Memc[image2], READ_ONLY, 0)

		len1 = IM_LENHDRMEM(im1) * SZ_STRUCT
		len2 = IM_LENHDRMEM(im2) * SZ_STRUCT
		ndx1 = len1 / SZ_RECORD + 1
		ndx2 = len2 / SZ_RECORD + 1

		call malloc (index1, ndx1, TY_INT)
		call malloc (index2, ndx2, TY_INT)
		call malloc (buffer1, len1, TY_CHAR)
		call malloc (buffer2, len2, TY_CHAR)

		call getkeybuf (im1, igroup, negate, Memc[keynames], 
				Memc[buffer1], Memi[index1], ndx1)
		call getkeybuf (im2, igroup, negate, Memc[keynames], 
				Memc[buffer2], Memi[index2], ndx2)

		call pr_hdiff (im1, im2, Memc[buffer1], Memi[index1], ndx1,
			       Memc[buffer2], Memi[index2], ndx2)

		call mfree (index1, TY_INT)
		call mfree (index2, TY_INT)
		call mfree (buffer1, TY_CHAR)
		call mfree (buffer2, TY_CHAR)

		call imunmap (im1)
		call imunmap (im2)
	    }
	    call tp_close (tp1)
	    call tp_close (tp2)
	}	
	call sfree (sp)
end
