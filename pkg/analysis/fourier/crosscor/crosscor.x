include	<imhdr.h>
include	"../fourier.h"
include	"../fterr.h"

# crosscor -- cross correlation
#
# Phil Hodge,  3-Aug-1988  Task created.
# Phil Hodge,  7-Feb-1990  Remove "tmp$" from names of temporary images.
# Phil Hodge,  4-Sep-1991  Normalize by dividing by the total number of pixels.
# Phil Hodge, 24-Jun-1992  Include verbose parameter.
# Phil Hodge,  2-Jul-1993  Do transforms of small 2-D arrays in-memory.
# Phil Hodge, 23-Jul-1993  Simplify computation of crpix, and move closer to
#			beginning of subroutine, before call to ft_small_x;
#			get in_memory before opening images.
# Phil Hodge, 15-Jan-1995  Call ft_struct_open and ft_struct_close.

procedure t_crosscor()
#--
pointer sp
pointer fti1, fti2, fto		# FT structures for input, output files
pointer infile1, infile2, outfile
int	len_blk			# edge length of block for transposing image
bool	getreal1, getreal2, getimag1, getimag2, crereal, creimag
bool	coord_shift		# include shift according to coordinates?
bool	center			# center output?
bool	chop			# truncate output to size of first input image?
bool	pad			# pad working area?
bool	in_memory		# do 2-D transforms in-memory?
bool	verbose			# print file names?
int	clgeti()
bool	clgetb()

begin
	call smark (sp)

	# Allocate Fourier structure for input & output images.
	call ft_struct_open (fti1)
	call ft_struct_open (fti2)
	call ft_struct_open (fto)

	call salloc (infile1, SZ_FNAME, TY_CHAR)
	call salloc (infile2, SZ_FNAME, TY_CHAR)
	call salloc (outfile, SZ_FNAME, TY_CHAR)

	# Get file names.
	call clgstr ("input1", Memc[infile1], SZ_FNAME)
	call clgstr ("input2", Memc[infile2], SZ_FNAME)
	call clgstr ("output", Memc[outfile], SZ_FNAME)

	if (Memc[outfile] == EOS)
	    call ft_error (FT_NOOUTPUT, FT_FATAL)

	# Get processing options.
	getreal1 = clgetb ("inreal1")
	getimag1 = clgetb ("inimag1")
	getreal2 = clgetb ("inreal2")
	getimag2 = clgetb ("inimag2")
	crereal = clgetb ("outreal")
	creimag = clgetb ("outimag")

	coord_shift = clgetb ("coord_shift")
	center = clgetb ("center")
	chop = clgetb ("chop")
	pad = clgetb ("pad")
	verbose = clgetb ("verbose")
	len_blk = 512				# default

	# Open images as requested.
	call ft_open_files_i (fti1, Memc[infile1], getreal1, getimag1)
	call ft_open_files_i (fti2, Memc[infile2], getreal2, getimag2)

	# Get in_memory if image is 2-D.
	if (FT_NAXIS(fti1) > 1) {
	    in_memory = clgetb ("inmemory")
	    if (in_memory && coord_shift)
		call error (1, "can't set inmemory=yes and coord_shift=yes")
	}

	# Create output image(s).  The first input image is used as the
	# template for creating the output image.
	call ft_open_files_o (fti1, fto, Memc[outfile], crereal, creimag)

	# Print file names.
	if (verbose) {
	    call ft_ri_print (fti1, "input1")
	    call ft_ri_print (fti2, "input2")
	    call ft_ri_print (fto, "output")
	    call flush (STDOUT)
	}

	# Perform cross correlation.
	switch (FT_NAXIS(fti1)) {
	case 1:
	    call ft_1_cross (fti1, fti2, fto, coord_shift, center, chop, pad)
	case 2:
	    len_blk = clgeti ("len_blk")
	    call ft_2_cross (fti1, fti2, fto,
			in_memory, len_blk, coord_shift, center, chop, pad)
	default:
	    call ft_error (FT_BADNAXIS, FT_FATAL)
	}

	if (FT_REPT(fto) != NULL) {
	    call sprintf (IM_TITLE(FT_REPT(fto)), SZ_IMTITLE,
			"Re (crosscor:  %s %s)")
		call pargstr (Memc[infile1])
		call pargstr (Memc[infile2])
	    call save_old_info (fti1, FT_REPT(fto))	# save crval, crpix
	    call save_ctstruct (fto, FT_REPT(fto))
	}

	if (FT_IMPT(fto) != NULL) {
	    call sprintf (IM_TITLE(FT_IMPT(fto)), SZ_IMTITLE,
			"Im (crosscor:  %s %s)")
		call pargstr (Memc[infile1])
		call pargstr (Memc[infile2])
	    call save_old_info (fti1, FT_IMPT(fto))
	    call save_ctstruct (fto, FT_IMPT(fto))
	}

	# Close files and clean up.
	call ft_close_files (fti1)
	call ft_close_files (fti2)
	call ft_close_files (fto)
	call sfree (sp)
	call ft_struct_close (fto)
	call ft_struct_close (fti1)
	call ft_struct_close (fti2)
end

# ft_1_cross -- 1-dimensional cross correlation

procedure ft_1_cross (fti1, fti2, fto, coord_shift, center, chop, pad)

pointer fti1		# i: FFT pointer; first input image
pointer fti2		# i: FFT pointer; second input image
pointer fto		# i: FFT pointer; output image
bool	coord_shift	# i: include shift according to coordinates?
bool	center		# i: center output?
bool	chop		# i: truncate output to size of first input image?
bool	pad		# i: pad working area?
#--
pointer work			# structure for work space and flags
pointer iRe1, iIm1		# ptrs for input data, first image (not hdr)
pointer iRe2, iIm2		# ptrs for input data, second image (not hdr)
pointer oRe, oIm		# ptrs to output data (not header struct)
pointer oRscr1, oIscr1		# for transform of 1st array; size = allpts
pointer oRscr2, oIscr2		# for output and transform of 2nd array
double	ref_crval[2]		# reference crval
double	ref_crpix[2]		# reference crpix
double	pshift[2]		# pixel offset between input images
real	xnpts			# =real(allpts)
int	npts1, npts2, allpts	# size of 1st & 2nd input, sum of both
int	outpts			# size of output (=npts1 or npts1+npts2)
int	cax			# current axis (one)
int	lowlim			# lower limit for copying to output
bool	fwd			# forward transform?  (first forward, then inv)
bool	not_zero[2]		# is phase shift non-trivial?

begin
	call malloc (work, LEN_WORK, TY_STRUCT)

	N_FILES(work) = 2		# two input file pairs
	fwd = true			# for sign of phase shift

	npts1 = IM_LEN(FT_IMAGE(fti1),1)
	npts2 = IM_LEN(FT_IMAGE(fti2),1)
	if (pad) 
		allpts = npts1 + npts2
	else
		allpts = max(npts1,npts2)
	NPTS(work) = allpts
	xnpts = real (allpts)

	# The size of the output image(s) depends on the CHOP parameter.
	if (chop)
	    outpts = npts1
	else
	    outpts = allpts

	if (FT_REPT(fto) != NULL)
	    IM_LEN(FT_REPT(fto),1) = outpts
	if (FT_IMPT(fto) != NULL)
	    IM_LEN(FT_IMPT(fto),1) = outpts

	# Copy coordinate info from input to output.  Also set FT_CRPIX = 1.
	call ft_wcs_n (fti1, fto)

	# Modify FT_CRPIX if output is to be centered.
	if (center)
	    FT_CRPIX(fto,1) = 1. + outpts / 2

	ref_crval[1] = FT_CRVAL(fti2,1)
	ref_crpix[1] = FT_CRPIX(fti2,1)
	call ft_phase_init_x (fti1, ref_crval, ref_crpix,
			pshift, not_zero)

	# Create scratch and/or read input.
	call ft_malloc (fti1, npts1, YES, iRe1, iIm1)
	call ft_malloc (fti2, npts2, YES, iRe2, iIm2)
	call ft_gl1r (FT_REPT(fti1), FT_IMPT(fti1), iRe1, iIm1)
	call ft_gl1r (FT_REPT(fti2), FT_IMPT(fti2), iRe2, iIm2)

	# Create scratch for transforms for input images.  The product
	# will be put in the second set of scratch arrays and later copied
	# to the output images.
	call malloc (oRscr1, allpts, TY_REAL)
	call malloc (oIscr1, allpts, TY_REAL)
	call malloc (oRscr2, allpts, TY_REAL)
	call malloc (oIscr2, allpts, TY_REAL)

	# For 1-D data, current axis is always 1.
	cax = 1

	# Initialization section.  The order is important here:
	# set SHIFT, call FT_WORK_INIT, call FT_PHASE_TAB.

	if (center)
	    SHIFT(work) = allpts / 2	# shift pixel 1 right by this amount

	# Create scratch for NCAR FFT routine, for phase shifting, centering.
	call ft_work_init (work, coord_shift, center)

	if (coord_shift)
	    call ft_phase_tab (work, pshift[1], not_zero[1])

	# Forward transform both input sets.
	call ft_fft_2 (work, npts1, npts2,
		Memr[iRe1], Memr[iIm1], Memr[iRe2], Memr[iIm2],
		Memr[oRscr1], Memr[oIscr1], Memr[oRscr2], Memr[oIscr2])

	# Phase shift?  It's applied to the transform of the first array.
	if (coord_shift && not_zero[1])
	    call ft_phase_s (allpts, Memr[oRscr1], Memr[oIscr1],
			Memr[COSTAB(work)], Memr[SINTAB(work)])

	# Cross correlation.  The output is put in scratch arrays #2.
	call ft_xcorr (Memr[oRscr1], Memr[oIscr1],
			Memr[oRscr2], Memr[oIscr2], allpts)

	# We're done with these scratch arrays.
	call mfree (oRscr1, TY_REAL)
	call mfree (oIscr1, TY_REAL)

	# Inverse transform of product (done in-place)
	fwd = false
	call ft_fft (work, Memr[oRscr2], Memr[oIscr2],
			   Memr[oRscr2], Memr[oIscr2], fwd)

	if (center) {
	    call ft_center (allpts, SHIFT(work),
			   Memr[oRscr2], Memr[oIscr2], Memr[C_COPY(work)])
	    if (chop)
		lowlim = 1 + allpts / 2 - outpts / 2
	    else
		lowlim = 1
	} else {
	    lowlim = 1
	}

	# Deallocate scratch for NCAR, etc (but not work itself).
	call ft_work_free (work, coord_shift, center)

	# Get pointers to output image(s).
	call ft_pl1r (FT_REPT(fto), FT_IMPT(fto), oRe, oIm)

	# Copy from scratch to output, normalize, and deallocate scratch.
	if (FT_REPT(fto) != NULL) {
	    call amovr (Memr[oRscr2+lowlim-1], Memr[oRe], outpts)
	    call adivkr (Memr[oRe], xnpts, Memr[oRe], outpts)
	    call mfree (oRscr2, TY_REAL)
	}
	if (FT_IMPT(fto) != NULL) {
	    call amovr (Memr[oIscr2+lowlim-1], Memr[oIm], outpts)
	    call adivkr (Memr[oIm], xnpts, Memr[oIm], outpts)
	    call mfree (oIscr2, TY_REAL)
	}

	# Deallocate memory buffers if we needed any.
	call ft_mfree (fti1, iRe1, iIm1)
	call ft_mfree (fti2, iRe2, iIm2)

	call mfree (work, TY_STRUCT)
end

# ft_2_cross -- 2-dimensional cross correlation
# Eight scratch images are used.  Four are for the transforms of the lines
# of the input file pairs (real & imaginary parts); four more are for the
# transpose of those.  This means 16 pointers:  eight for the imhdr
# struct and eight pointers to data.  The pointers to the imhdr struct
# all contain an "_", while those to data do not.  An "x" indicates the
# pointer refers to the transposed image.  An "s" prefix indicates a pointer
# to data for a scratch image.

procedure ft_2_cross (fti1, fti2, fto,
		in_memory, len_blk, coord_shift, center, chop, pad)

pointer fti1		# i: FFT pointer; input image
pointer fti2		# i: FFT pointer; input image
pointer fto		# i: FFT pointer; output image
bool	in_memory	# i: do transforms in-memory?
int	len_blk		# i: block size for transposing image
bool	coord_shift	# i: include shift according to coordinates?
bool	center		# i: center output?
bool	chop		# i: truncate output to size of first input image?
bool 	pad		# i: pad working area?
#--
pointer work			# structure for work space and flags
pointer iRe1, iIm1		# ptrs for input data, first image (not hdr)
pointer iRe2, iIm2		# ptrs for input data, second image (not hdr)
pointer oRe, oIm		# ptrs to output data (not header struct)
pointer Re_1, Im_1		# imhdr struct for transform of 1st image
pointer Re_2, Im_2		# imhdr struct for transform of 2nd image
pointer Re_x1, Im_x1		# imhdr struct for 1st transposed images
pointer Re_x2, Im_x2		# imhdr struct for 2nd transposed images
pointer sRe1, sIm1		# ptr to data for transform of 1st image
pointer sRe2, sIm2		# ptr to data for transform of 2nd image
pointer sRex1, sImx1		# ptr to data for 1st transposed image
pointer sRex2, sImx2		# ptr to data for 2nd transposed image
double	ref_crval[2]		# reference crval
double	ref_crpix[2]		# reference crpix
double	pshift[2]		# pixel offset between input images
real	xnpts2			# = allpts[1] * allpts[2], for normalizing
int	npts1[2], npts2[2]	# size of 1st & 2nd input
int	n1, n2			# npts1[1] & npts2[1] or zero if outside image
int	allpts[2]		# npts1 + npts2 = size of scratch & transposed
int	outpts[2]		# size of output (=npts1 or npts1+npts2)
int	lowlim[2]		# lower limit for copying to output
int	cax			# current axis
int	k			# loop index
int	line, column		# loop indexes for line & column of image
bool	fwd			# forward transform?
bool	not_zero[2]		# is phase shift non-trivial?
string	seed1r "ft1r"		# for transform of first input
string	seed1i "ft1i"
string	seed2r "ft2r"		# for transform of second input
string	seed2i "ft2i"
string	sedx1r "ftx1r"		# for transpose of first
string	sedx1i "ftx1i"
string	sedx2r "ftx2r"		# for transpose of second
string	sedx2i "ftx2i"

begin
	# Copy coordinate info from input to output.
	call ft_wcs_n (fti1, fto)

	# Get offset between the two input images based on crval, etc.
	do k = 1, 2 {
	    ref_crval[k] = FT_CRVAL(fti2,k)
	    ref_crpix[k] = FT_CRPIX(fti2,k)
	}
	call ft_phase_init_x (fti1, ref_crval, ref_crpix, pshift, not_zero)

	fwd = true

	# For each axis, get size of input, scratch & transposed, output files.
	do k = 1, 2 {
	    npts1[k] = IM_LEN(FT_IMAGE(fti1),k)
	    npts2[k] = IM_LEN(FT_IMAGE(fti2),k)
	    if (pad)
		allpts[k] = npts1[k] + npts2[k]
	    else
		allpts[k] = max(npts1[k],npts2[k])
	    if (chop)
		outpts[k] = npts1[k]
	    else
		outpts[k] = allpts[k]
	    if (FT_REPT(fto) != NULL)
		IM_LEN(FT_REPT(fto),k) = outpts[k]
	    if (FT_IMPT(fto) != NULL)
		IM_LEN(FT_IMPT(fto),k) = outpts[k]
	}
	xnpts2 = real (allpts[1]) * real (allpts[2])	# for normalizing

	# Modify FT_CRPIX if output is to be centered.
	# The call to ft_wcs_n has already set FT_CRPIX = 1.
	if (center) {
	    FT_CRPIX(fto,1) = 1 + outpts[1] / 2		# truncate
	    FT_CRPIX(fto,2) = 1 + outpts[2] / 2
	}

	# Do the transform in memory?
	if (in_memory) {
	    call ft_small_x (fti1, fti2, fto,
			allpts[1], allpts[2], center)
	    return
	}

	call malloc (work, LEN_WORK, TY_STRUCT)
	N_FILES(work) = 2		# two input file pairs

	# Create scratch files for the transforms of the lines of the
	# two input files.
	call ft_map_s (seed1r, 2, allpts[1], allpts[2], Re_1)
	call ft_map_s (seed1i, 2, allpts[1], allpts[2], Im_1)
	call ft_map_s (seed2r, 2, allpts[1], allpts[2], Re_2)
	call ft_map_s (seed2i, 2, allpts[1], allpts[2], Im_2)

	# ********************
	# Begin the section for taking the transform of each line.

	cax = 1					# first axis
	NPTS(work) = allpts[cax]

	# Allocate scratch for input files that don't exist.
	call ft_malloc (fti1, npts1[cax], YES, iRe1, iIm1)
	call ft_malloc (fti2, npts2[cax], YES, iRe2, iIm2)

	# Initialization section.
	# Create scratch for NCAR FFT routine and for phase shifting.
	# We cannot center the output yet.
	call ft_work_init (work, coord_shift, NO)

	# Compute cosines & sines for first axis.
	if (coord_shift)
	    call ft_phase_tab (work, pshift[1], not_zero[1])

	# Forward transform each line of each input set.
	n1 = npts1[cax]
	n2 = npts2[cax]
	do line = 1, allpts[2] {		# do for each line
	    if (line <= npts1[2])
		call ft_gl2r (FT_REPT(fti1), FT_IMPT(fti1), line, iRe1, iIm1)
	    else
		n1 = 0				# indicates no data

	    if (line <= npts2[2])
		call ft_gl2r (FT_REPT(fti2), FT_IMPT(fti2), line, iRe2, iIm2)
	    else
		n2 = 0				# no data

	    call ft_pl2r (Re_1, Im_1, line, sRe1, sIm1)
	    call ft_pl2r (Re_2, Im_2, line, sRe2, sIm2)

	    call ft_fft_2 (work, n1, n2,
		Memr[iRe1], Memr[iIm1], Memr[iRe2], Memr[iIm2],
		Memr[sRe1], Memr[sIm1], Memr[sRe2], Memr[sIm2])

	    # Phase shift?  It's applied to the transform of the first array.
	    if (coord_shift  &&  not_zero[1]  &&  n1>0)
		call ft_phase_s (allpts[cax], Memr[sRe1], Memr[sIm1],
			Memr[COSTAB(work)], Memr[SINTAB(work)])
	}

	# Deallocate scratch for NCAR, etc (but not work itself).
	call ft_work_free (work, coord_shift, NO)

	# Flush scratch images so we can read them.
	call imflush (Re_1)
	call imflush (Im_1)
	call imflush (Re_2)
	call imflush (Im_2)

	# Deallocate memory buffers for input if we needed any.
	call ft_mfree (fti1, iRe1, iIm1)
	call ft_mfree (fti2, iRe2, iIm2)

	# ********************
	# Create scratch files for the transpose of the scratch image.
	# (Note that the lengths of the axes are interchanged.)
	# Transpose the scratch file, and delete the first scratch file.
	# We still need Re_2 & Im_2.

	call ft_map_s (sedx1r, 2, allpts[2], allpts[1], Re_x1)
	call ft_transpose (Re_1, Re_x1, len_blk)
	call ft_del_s (Re_1)

	call ft_map_s (sedx1i, 2, allpts[2], allpts[1], Im_x1)
	call ft_transpose (Im_1, Im_x1, len_blk)
	call ft_del_s (Im_1)

	call ft_map_s (sedx2r, 2, allpts[2], allpts[1], Re_x2)
	call ft_map_s (sedx2i, 2, allpts[2], allpts[1], Im_x2)

	call ft_transpose (Re_2, Re_x2, len_blk)
	call ft_transpose (Im_2, Im_x2, len_blk)

	# Flush so we can read & write.
	call imflush (Re_x1)
	call imflush (Im_x1)
	call imflush (Re_x2)
	call imflush (Im_x2)

	# ********************
	# Begin the section for taking the transform of each column.

	cax = 2					# second axis
	NPTS(work) = allpts[cax]

	# Initialization section.
	if (center)
	    SHIFT(work) = allpts[cax] / 2

	call ft_work_init (work, coord_shift, center)

	# Compute cosines & sines for second axis.
	if (coord_shift)
	    call ft_phase_tab (work, pshift[2], not_zero[2])

	# Forward transform both sets in-place in the transposed images,
	# multiply in order to do the cross correlation, and then inverse
	# transform the columns.  (We still have to inv. transform the rows.)

	do column = 1, allpts[1] {		# do for each column

	    call ft_gl2r (Re_x1, Im_x1, column, sRex1, sImx1)
	    call ft_gl2r (Re_x2, Im_x2, column, sRex2, sImx2)
	    call ft_pl2r (Re_x1, Im_x1, column, sRe1, sIm1)
	    call ft_pl2r (Re_x2, Im_x2, column, sRe2, sIm2)

	    # This is a forward transform.
	    call ft_fft_2 (work, allpts[cax], allpts[cax],
		Memr[sRex1], Memr[sImx1], Memr[sRex2], Memr[sImx2],
		Memr[sRe1], Memr[sIm1], Memr[sRe2], Memr[sIm2])

	    # Phase shift?  It's applied to the transform of the first array.
	    if (coord_shift && not_zero[2])
		call ft_phase_s (allpts[cax], Memr[sRe1], Memr[sIm1],
			Memr[COSTAB(work)], Memr[SINTAB(work)])

	    # Cross correlation.  Scratch arrays #2 are input & output.
	    call ft_xcorr (Memr[sRe1], Memr[sIm1],
			Memr[sRe2], Memr[sIm2], allpts[cax])

	    # Inverse transform of product (done in-place).
	    fwd = false
	    call ft_fft (work, Memr[sRe2], Memr[sIm2],
			   Memr[sRe2], Memr[sIm2], fwd)

	    if (center)
		call ft_center (allpts[cax], SHIFT(work),
			   Memr[sRe2], Memr[sIm2], Memr[C_COPY(work)])
	}
	call ft_work_free (work, coord_shift, center)

	# ********************
	# Transpose back to the scratch files, and delete the scratch
	# transposed files.
	call ft_del_s (Re_x1)			# don't need #1 any more
	call ft_del_s (Im_x1)

	call imflush (Re_x2)
	call imflush (Im_x2)

	call ft_transpose (Re_x2, Re_2, len_blk)
	call ft_del_s (Re_x2)

	call ft_transpose (Im_x2, Im_2, len_blk)
	call ft_del_s (Im_x2)

	call imflush (Re_2)
	call imflush (Im_2)

	# ********************
	# Now return to the first axis and take the inverse transform.
	# In this section we also copy the result to the output file(s).

	cax = 1
	NPTS(work) = allpts[cax]

	# Initialization section.
	if (center)
	    SHIFT(work) = allpts[cax] / 2

	# Get the lower limits for extracting the portion of the data
	# that we want to keep, as specified by CENTER and CHOP.
	do k = 1, 2 {
	    if (center && chop)
		lowlim[k] = 1 + allpts[k] / 2 - outpts[k] / 2
	    else
		lowlim[k] = 1
	}

	# No coordinate shifting remains to be done.
	call ft_work_init (work, NO, center)

	fwd = false				# inverse transform
	# Do for each line that is to be written to output.
	line = 1
	do k = lowlim[2], lowlim[2]+outpts[2]-1 {

	    call ft_gl2r (Re_2, Im_2, k, sRe2, sIm2)		# get line K

	    call ft_fft (work, Memr[sRe2], Memr[sIm2],
			   Memr[sRe2], Memr[sIm2], fwd)		# in-place

	    if (center)
		call ft_center (allpts[cax], SHIFT(work),
			   Memr[sRe2], Memr[sIm2], Memr[C_COPY(work)])

	    # Copy to output; put at line number LINE; normalize.
	    call ft_pl2r (FT_REPT(fto), FT_IMPT(fto), line, oRe, oIm)
	    if (FT_REPT(fto) != NULL) {
		call amovr (Memr[sRe2+lowlim[1]-1], Memr[oRe], outpts[1])
		call adivkr (Memr[oRe], xnpts2, Memr[oRe], outpts[1])
	    }
	    if (FT_IMPT(fto) != NULL) {
		call amovr (Memr[sIm2+lowlim[1]-1], Memr[oIm], outpts[1])
		call adivkr (Memr[oIm], xnpts2, Memr[oIm], outpts[1])
	    }

	    line = line + 1
	}

	call ft_work_free (work, NO, center)

	call ft_del_s (Re_2)
	call ft_del_s (Im_2)

	call mfree (work, TY_STRUCT)
end
