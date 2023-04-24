include	<error.h>
include	"addmasks.h"
define	done_		99

# ADDMASKS -- Combine several masks
#
# B.Simon	04-Apr-91	Original
# B.Simon	16-Jul-93	Revised to use expression evaluator

procedure addmasks ()

#--
pointer	input		# list of input mask names
pointer	output		# output mask name
pointer	expr		# expression to evaluate
pointer	flags		# list of flag values

int	ifile, nfile, nflag, junk
int	iline, len1, len2, nline1, nline2
pointer	sp, errmsg, file, fname, flag, order, pcode
pointer	inmask, outmask, inbuf, outbuf

string	nomasks  "Input mask not found (%s)"
string	badsize  "Mask size does not match first mask (%s)"

int	imtlen(), imtgetim()
pointer	imtopen(), addcompile()
pointer	old_mask(), new_mask(), rd_mask(), wrt_mask()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (expr, SZ_FNAME, TY_CHAR)
	call salloc (flags, SZ_FNAME, TY_CHAR)

	call salloc (fname, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)
	call salloc (flag, MAX_FLAGS, TY_INT)
	call salloc (order, MAX_FLAGS, TY_INT)

	# Read task parameters

	call clgstr ("input", Memc[input], SZ_FNAME)
	call clgstr ("output", Memc[output], SZ_FNAME)
	call clgstr ("expr", Memc[expr], SZ_FNAME)
	call clgstr ("flags", Memc[flags], SZ_FNAME)

	# Compile the expression

	pcode = addcompile (Memc[expr])

	# Decode list of flag values

	call rd_flags (Memc[flags], Memi[flag], nflag, MAX_FLAGS)
	call sort_flag (Memi[flag], Memi[order], nflag)

	# Open input mask template

	file = imtopen (Memc[input])
	nfile = imtlen (file)
	if (nfile == 0) {
	    call sprintf (Memc[errmsg], SZ_LINE, nomasks)
	    call pargstr (Memc[input])
	    call error (1, Memc[errmsg])
	}

	# Initialize warning message subroutine with filenames

	call set_warn (file, output)

	# Open input masks. Check sizes to see that they all agree

	call salloc (inmask, nfile, TY_INT)
	call salloc (inbuf, nfile, TY_INT)

	do ifile = 1, nfile {
	    junk = imtgetim (file, Memc[fname], SZ_FNAME)
	    Memi[inmask+ifile-1] = old_mask (Memc[fname], READ_ONLY)

	    if (ifile == 1) {
		call sz_mask (Memi[inmask], len1, nline1)
	    } else {
		call sz_mask (Memi[inmask+ifile-1], len2, nline2)
		if ((len1 != len2) || (nline1 != nline2)) {
		    call sprintf (Memc[errmsg], SZ_LINE, badsize)
		    call pargstr (Memc[fname])
		    call error (1, Memc[errmsg])
		}
	    }
	}

	# Open output mask

	outmask = new_mask (Memc[output], 0, Memi[inmask])

	# Loop over each line in the mask

	iline = 0
	repeat {
	    iline = iline + 1

	    # Read next line from input files

	    do ifile = 1, nfile {
		Memi[inbuf+ifile-1] = rd_mask (Memi[inmask+ifile-1])
		if (Memi[inbuf+ifile-1] == NULL)
		    goto done_
	    }

	    outbuf = wrt_mask (outmask)

	    # Evaluate expression and store result in output buffer

	    iferr {
		call addeval (pcode, iline, Memi[flag], Memi[order], 
			      nflag, Memi[inbuf], outbuf, len1, nfile)

	    } then {
		call del_mask (outmask)
		call erract (EA_ERROR)
	    }
	}

	# Close mask files

done_	do ifile = 1, nfile
	    call cls_mask (Memi[inmask+ifile-1])

	call cls_mask (outmask)
	call imtclose (file)

	# Free dynamic memory

	call addfree (pcode)
	call sfree (sp)
end
