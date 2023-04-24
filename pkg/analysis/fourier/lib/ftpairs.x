include	<imhdr.h>
include "../fourier.h"
define	SZ_CTYPE2	(SZ_CTYPE + 2)	# space for trailing blanks

# FT_READ_FTPAIRS -- Read ASCII list of CTYPE transform pairs.
#
# CTYPE is transformed via a lookup table.  The LUT data are stored in an
# ASCII file that is read into memory by this routine.  There is a set of
# default CTYPE transformations that is based on the "recommended" CTYPE names
# in the FITS papers.  The name of the file is a task parameter, and its
# format is pretty obvious, so users may create their own file if they have
# data in a coordinate system that is not covered by the default list.
#
# Memory is allocated for the CTYPE dictionary, and the pointers are
# returned in the two-element array CTD.  CTD can be deallocated by calling
# ft_pairs_free.
#
# C. D. Biemesderfer, STScI, 21 Dec 87
# Phil Hodge, 30-Dec-1988  Realloc if dictionary is full.

procedure ft_read_ftpairs (ctd)

pointer ctd[2]			# o: pointers to input & output CTYPE arrays
#--
pointer sp
pointer ftp			# scratch for name of ftpairs file
char	ctype1[SZ_CTYPE2], ctype2[SZ_CTYPE2], blank[SZ_CTYPE2]
int	sz_dic			# default size of CTYPE dictionary
int	k			# loop index
int	np			# number of pairs
int	ftpd			# fd for ftpairs file
int	open(), fscan()

begin
	call smark (sp)
	call salloc (ftp, SZ_FNAME, TY_CHAR)

	do k = 1, SZ_CTYPE2
	    blank[k] = ' '
	blank[SZ_CTYPE2+1] = EOS

	# Allocate space for the ctype dictionary.
	sz_dic = MAX_CTPAIRS * SZ_CTYPE2	# initial size
	call malloc (ctd[1], sz_dic, TY_CHAR)
	call malloc (ctd[2], sz_dic, TY_CHAR)

	# The first CTYPE pair is always the null pair.
	Memc[ctd[1]] = EOS
	Memc[ctd[2]] = EOS
	call strcat (blank, Memc[ctd[1]], sz_dic)
	call strcat (blank, Memc[ctd[2]], sz_dic)
	np = 1					# one pair so far

	# Open text file containing CTYPE pairs.
	call clgstr ("ftpairs", Memc[ftp], SZ_FNAME)
	iferr (ftpd = open (Memc[ftp], READ_ONLY, TEXT_FILE)) {
	    call eprintf ("FT pkg:  can't open ftpairs file:\n   %s\n")
		call pargstr (Memc[ftp])
	    call sfree (sp)
	    return
	}

	# Read through ftpairs file and load CTYPE dictionaries.  Null records
	# and comments are skipped.  Original and transform domain CTYPEs are
	# expected as the first two white-space delimited fields on each record.

	while (fscan (ftpd) != EOF) {

	    call gargwrd (ctype1, SZ_CTYPE)	# SZ_CTYPE, not SZ_CTYPE2
	    call gargwrd (ctype2, SZ_CTYPE)
	    if (ctype1[1] == EOS || ctype1[1] == '#')
		next

	    np = np + 1
	    if (np >= sz_dic) {
		# Double the size of the input & output dictionaries.
		sz_dic = 2 * sz_dic
		call realloc (ctd[1], sz_dic, TY_CHAR)
		call realloc (ctd[2], sz_dic, TY_CHAR)
	    }

	    call strcat (blank, ctype1, SZ_CTYPE2)
	    call strcat (ctype1, Memc[ctd[1]], sz_dic)	# input CTYPE

	    if (ctype2[1] == EOS || ctype2[1] == '#')
		call strcpy (blank, ctype2, SZ_CTYPE2)
	    else
		call strcat (blank, ctype2, SZ_CTYPE2)
	    call strcat (ctype2, Memc[ctd[2]], sz_dic)	# output CTYPE
	}

	call close (ftpd)
	call sfree (sp)
end

# FT_CHANGE_CTYPE -- Determine output (transformed) coordinate type.
#
# The input dictionary is searched for the given CTYPE, and the transform
# domain CTYPE set to the corresponding entry in the output dictionary.

procedure ft_change_ctype (fwd, ctdin, ctdout, inct, outct)

bool	fwd		# i: forward transform?
char	ctdin[ARB]	# i: buffer for input CTYPEs
char	ctdout[ARB]	# i: buffer for output CTYPEs
char	inct[ARB]	# i: input CTYPE
char	outct[ARB]	# o: CTYPE in transform domain
#--
int	octbeg, octend
int	gstrmatch()

begin
	if ( fwd ) {
	    if (gstrmatch (ctdin, inct, octbeg, octend) > 0)
		call strcpy (ctdout[octbeg], outct, SZ_CTYPE)
	    else
		call strcpy ("PIXEL", outct, SZ_CTYPE)
	} else {
	    if (gstrmatch (ctdout, inct, octbeg, octend) > 0)
		call strcpy (ctdin[octbeg], outct, SZ_CTYPE)
	    else
		call strcpy ("PIXEL", outct, SZ_CTYPE)
	}
end

# FT_PAIRS_FREE -- free memory for CTYPE dictionary
# This routine releases memory that was allocated for the lists of
# input & output CTYPEs.

procedure ft_pairs_free (ctd)

pointer ctd[2]		# io: pointers to memory for CTYPE lists
#--

begin
	if (ctd[2] != NULL)
	    call mfree (ctd[2], TY_CHAR)
	if (ctd[1] != NULL)
	    call mfree (ctd[1], TY_CHAR)
end
