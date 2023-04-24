include	<imhdr.h>
include	"../fourier.h"
include	"../fterr.h"

# inverse -- inverse Fourier transform
#
# C. D. Biemesderfer, STScI, 21 Dec 87
# Phil Hodge, 16-Aug-1988  Add 2-D.
# Phil Hodge,  7-Feb-1990  Remove "tmp$" from names of temporary images.
# Phil Hodge, 19-Apr-1992  Call ft_small_i for 2-D small images; incl verbose.
# Phil Hodge, 23-Jul-1993  Get in_memory before opening images.
# Phil Hodge, 15-Jan-1996  Call ft_struct_open and ft_struct_close.

procedure t_inverse()
#--
pointer sp
pointer fti, fto		# FT structures for input, output files
pointer ctd[2]			# pointers to input & output CTYPE dictionaries
pointer infile, outfile
int	len_blk			# edge length of block for transposing image
bool	getreal, getimag, crereal, creimag
bool	coord_shift		# include shift according to coordinates?
bool	decenter		# decenter input?
bool	in_memory		# do 2-D transform in-memory?
bool	verbose			# print file names?
int	clgeti()
bool	clgetb()

begin
	call smark (sp)

	# Allocate Fourier structure for input & output images.
	call ft_struct_open (fti)
	call ft_struct_open (fto)

	call salloc (infile, SZ_FNAME, TY_CHAR)
	call salloc (outfile, SZ_FNAME, TY_CHAR)

	# Get file names.
	call clgstr ("input", Memc[infile], SZ_FNAME)
	call clgstr ("output", Memc[outfile], SZ_FNAME)

	if (Memc[outfile] == EOS)
	    call ft_error (FT_NOOUTPUT, FT_FATAL)

	# Get processing options.
	getreal = clgetb ("inreal")
	getimag = clgetb ("inimag")
	crereal = clgetb ("outreal")
	creimag = clgetb ("outimag")

	coord_shift = clgetb ("coord_shift")
	decenter = clgetb ("decenter")		# enable decentering of input?
	in_memory = clgetb ("inmemory")
	verbose = clgetb ("verbose")

	# Get "ftpairs" parameter, and read that file.
	call ft_read_ftpairs (ctd)

	# Open input image(s).
	call ft_open_files_i (fti, Memc[infile], getreal, getimag)

	# Get in_memory if image is 2-D.
	if (FT_NAXIS(fti) > 1) {
	    in_memory = clgetb ("inmemory")
	    if (in_memory && coord_shift)
		call error (1, "can't set inmemory=yes and coord_shift=yes")
	}

	# Create output image(s).
	call ft_open_files_o (fti, fto, Memc[outfile], crereal, creimag)

	# Print file names.
	if (verbose) {
	    call ft_ri_print (fti, "input")
	    call ft_ri_print (fto, "output")
	    call flush (STDOUT)
	}

	# Perform transform.
	switch (FT_NAXIS(fti)) {
	case 1:
	    call ft_1d_inverse (fti, fto, ctd, coord_shift, decenter)
	case 2:
	    len_blk = clgeti ("len_blk")
	    call ft_2d_inverse (fti, fto, ctd, coord_shift, decenter,
			in_memory, len_blk)
	default:
	    call ft_error (FT_BADNAXIS, FT_FATAL)
	}

	if (FT_REPT(fto) != NULL) {
	    call sprintf (IM_TITLE(FT_REPT(fto)), SZ_IMTITLE,
			"Re (inv FT of %s)")
		call pargstr (Memc[infile])
	    call save_old_info (fti, FT_REPT(fto))	# save crval, crpix
	    call save_ctstruct (fto, FT_REPT(fto))
	}

	if (FT_IMPT(fto) != NULL) {
	    call sprintf (IM_TITLE(FT_IMPT(fto)), SZ_IMTITLE,
			"Im (inv FT of %s)")
		call pargstr (Memc[infile])
	    call save_old_info (fti, FT_IMPT(fto))
	    call save_ctstruct (fto, FT_IMPT(fto))
	}

	# Close files and clean up.
	call ft_close_files (fti)
	call ft_close_files (fto)
	call sfree (sp)
	call ft_pairs_free (ctd)
	call ft_struct_close (fto)
	call ft_struct_close (fti)
end

# ft_1d_inverse -- 1-dimensional inverse transform.

procedure ft_1d_inverse (fti, fto, ctd, coord_shift, decenter)

pointer fti		# i: ptr to FT struct for input image
pointer fto		# i: ptr to FT struct for output image
pointer ctd[2]		# i: pointers to input & output CTYPE dictionaries
bool	coord_shift	# i: include shift according to coordinates?
bool	decenter	# i: decenter input?
#--
pointer work		# structure for work space and flags
pointer iRe, iIm
pointer oRe, oIm
double	pshift[2]		# offset of crpix from 1
real	xnpts			# = npts
int	npts
int	cax			# current axis (one)
bool	fwd			# forward transform?  (no)
bool	not_zero[2]		# is phase shift non-trivial?

begin
	call malloc (work, LEN_WORK, TY_STRUCT)

	fwd = false			# inverse transform
	N_FILES(work) = 1		# one input file pair

	# Assign coordinate info in fto.  coord_shift will be reset
	# to false if ocrval & ocrpix not found.
	call ft_wcs_inv (fti, fto, ctd, coord_shift)

	# Get the offset from the first pixel to the reference pixel.
	call ft_phase_init (fto, pshift, not_zero)
	pshift[1] = -pshift[1]		# inverse transform

	npts = IM_LEN(FT_IMAGE(fti),1)
	NPTS(work) = npts

	# Get scratch space in place of files that don't exist.
	call ft_malloc (fti, npts, YES, iRe, iIm)	# yes, init to zero
	call ft_malloc (fto, npts,  NO, oRe, oIm)

	# Get pointers to files that do exist.
	call ft_gl1r (FT_REPT(fti), FT_IMPT(fti), iRe, iIm)
	call ft_pl1r (FT_REPT(fto), FT_IMPT(fto), oRe, oIm)

	cax = 1			# for 1-D data current axis is always 1

	# Initialization section.  The order is important here:
	# set SHIFT, call FT_WORK_INIT, call FT_PHASE_TAB.

	if (decenter)
	    # Shift pixel 1 (right or left, depending on sign) by this amount.
	    SHIFT(work) = 1 - nint (FT_CRPIX(fti,cax))

	# Create scratch for NCAR FFT routine, for phase shifting, decentering.
	call ft_work_init (work, coord_shift, decenter)
	if (coord_shift)
	    call ft_phase_tab (work, pshift[1], not_zero[1])

	if (decenter)
	    call ft_center (npts, SHIFT(work),
			Memr[iRe], Memr[iIm], Memr[C_COPY(work)])

	if (coord_shift && not_zero[1])
	    call ft_phase_s (npts, Memr[iRe], Memr[iIm],
			Memr[COSTAB(work)], Memr[SINTAB(work)])

	# Inverse FFT.
	call ft_fft (work, Memr[iRe], Memr[iIm], Memr[oRe], Memr[oIm], fwd)

	# Normalization.  The forward and inverse transforms and the power
	# spectrum are mutually normalized so that the answers among the
	# three of them are consistent.  The convention adopted is such that
	# the inverse transform requires a 1/N factor.

	xnpts = real(npts)
	call adivkr (Memr[oRe], xnpts, Memr[oRe], npts)
	call adivkr (Memr[oIm], xnpts, Memr[oIm], npts)

	# Deallocate scratch for NCAR, etc (but not work itself).
	call ft_work_free (work, coord_shift, decenter)

	# Deallocate memory buffers if we needed any.
	call ft_mfree (fto, oRe, oIm)
	call ft_mfree (fti, iRe, iIm)

	call mfree (work, TY_STRUCT)
end

# ft_2d_inverse -- 2-D inverse Fourier transform
# Both the real & imaginary parts for the output file must be created.
# The unwanted part (if either) will be deleted when we're done, however.

procedure ft_2d_inverse (fti, fto, ctd, coord_shift, decenter,
		in_memory, len_blk)

pointer fti		# i: ptr to FT struct for input image
pointer fto		# i: ptr to FT struct for output image
pointer ctd[2]		# i: pointers to input & output CTYPE dictionaries
bool	coord_shift	# i: include shift according to coordinates?
bool	decenter	# i: decenter input?
bool	in_memory	# i: do 2-D transform in-memory?
int	len_blk		# i: block size for transposing image
#--
pointer work		# structure for work space and flags
pointer iRe, iIm	# input
pointer oRe, oIm	# output
pointer Re_scr, Im_scr	# ptr to imhdr structure for scratch images
pointer ip		# ptr to imhdr structure for input image
double	pshift[2]	# offset of crpix from (1,1)
real	xnpts2		# = npts1 * npts2
int	npts		# length of current axis
int	cax		# current axis
int	line		# loop index for line of image
bool	crereal		# create output real part?
bool	creimag		# create output imaginary part?
bool	fwd		# forward transform?  (no)
bool	not_zero[2]	# is phase shift non-trivial?
pointer imgl2r(), impl2r()
string	seedri "ftri"
string	seedr  "ftr"
string	seedi  "fti"

begin
	# Assign coordinate info in fto.  coord_shift will be reset
	# to false if ocrval & ocrpix are not found.
	call ft_wcs_inv (fti, fto, ctd, coord_shift)

	ip = FT_IMAGE(fti)

	# Do the transform in memory?
	if (in_memory) {
	    call ft_small_i (fti, fto, decenter)
	    return
	}

	call malloc (work, LEN_WORK, TY_STRUCT)

	fwd = false			# inverse transform
	N_FILES(work) = 1		# one input file pair

	# Get the offset from the first pixel to the reference pixel.
	call ft_phase_init (fto, pshift, not_zero)
	pshift[1] = -pshift[1]
	pshift[2] = -pshift[2]

	# Both real & imaginary output files must exist while we're doing the
	# transform, but we can delete either afterwards.  Save the flags,
	# and create a scratch file if either part hasn't been created yet.
	crereal = (FT_REPT(fto) != NULL)
	creimag = (FT_IMPT(fto) != NULL)
	if ( ! crereal ) {
	    FT_REAL(fto) = YES			# reset later
	    call ft_map_s (seedri, 2, IM_LEN(ip,1), IM_LEN(ip,2), FT_REPT(fto))
	} else if ( ! creimag ) {
	    FT_IMAG(fto) = YES			# reset later
	    call ft_map_s (seedri, 2, IM_LEN(ip,1), IM_LEN(ip,2), FT_IMPT(fto))
	}

	# Begin the section for taking the transform of each line.
	# Write the results into the output file, but at this point
	# we're just using it for scratch.

	cax = 1					# first axis
	npts = IM_LEN(ip,cax)
	NPTS(work) = npts

	# Allocate scratch for input files that don't exist.
	call ft_malloc (fti, npts, YES, iRe, iIm)

	# Initialization section for work struct for an array of size npts.
	# The order is important here:
	# set SHIFT, call FT_WORK_INIT, call FT_PHASE_TAB.

	if (decenter)
	    SHIFT(work) = 1 - nint (FT_CRPIX(fti,cax))

	# Create scratch for NCAR FFT routine, for phase shifting, decentering.
	call ft_work_init (work, coord_shift, decenter)

	if (coord_shift)
	    call ft_phase_tab (work, pshift[1], not_zero[1])

	do line = 1, IM_LEN(ip,2) {		# do for each line

	    call ft_gl2r (FT_REPT(fti), FT_IMPT(fti), line, iRe, iIm)
	    oRe = impl2r (FT_REPT(fto), line)
	    oIm = impl2r (FT_IMPT(fto), line)

	    if (decenter)
		call ft_center (npts, SHIFT(work),
			Memr[iRe], Memr[iIm], Memr[C_COPY(work)])

	    if (coord_shift && not_zero[1])
		call ft_phase_s (npts, Memr[iRe], Memr[iIm],
			Memr[COSTAB(work)], Memr[SINTAB(work)])

	    # inverse FFT.
	    call ft_fft (work, Memr[iRe], Memr[iIm], Memr[oRe], Memr[oIm], fwd)
	}

	# Flush output so we can read from it.
	call imflush (FT_REPT(fto))
	call imflush (FT_IMPT(fto))

	# Deallocate scratch for NCAR, etc (but not work itself)
	call ft_work_free (work, coord_shift, decenter)

	# Deallocate memory buffers if we needed any.
	call ft_mfree (fti, iRe, iIm)

	# ********************
	# Create scratch files for the transpose of the output image.  The
	# lengths of the axes are interchanged; Re_scr & Im_scr are output.
	call ft_map_s (seedr, 2, IM_LEN(ip,2), IM_LEN(ip,1), Re_scr)
	call ft_map_s (seedi, 2, IM_LEN(ip,2), IM_LEN(ip,1), Im_scr)

	# Transpose the output file into the scratch file.
	call ft_transpose (FT_REPT(fto), Re_scr, len_blk)
	call ft_transpose (FT_IMPT(fto), Im_scr, len_blk)

	# Flush transposed image so we can read from it.
	call imflush (Re_scr)
	call imflush (Im_scr)

	# ********************
	# Begin the section for taking the transform of each column.

	cax = 2					# second axis
	npts = IM_LEN(ip,cax)
	NPTS(work) = npts

	# Initialization section for work struct.

	if (decenter)
	    SHIFT(work) = 1 - nint (FT_CRPIX(fti,cax))

	call ft_work_init (work, coord_shift, decenter)

	if (coord_shift)
	    call ft_phase_tab (work, pshift[2], not_zero[2])

	# We will divide by npts1*npts2 to normalize.
	xnpts2 = real(IM_LEN(ip,1)) * real(IM_LEN(ip,2)) 

	do line = 1, IM_LEN(ip,1) {		# do for each column

	    # We're doing this in-place in the scratch (transposed) images.
	    iRe = imgl2r (Re_scr, line)
	    iIm = imgl2r (Im_scr, line)
	    oRe = impl2r (Re_scr, line)
	    oIm = impl2r (Im_scr, line)

	    if (decenter)
		call ft_center (npts, SHIFT(work),
			Memr[iRe], Memr[iIm], Memr[C_COPY(work)])

	    if (coord_shift && not_zero[2])
		call ft_phase_s (npts, Memr[iRe], Memr[iIm],
			Memr[COSTAB(work)], Memr[SINTAB(work)])

	    call ft_fft (work, Memr[iRe], Memr[iIm], Memr[oRe], Memr[oIm], fwd)

	    # Normalize the output.
	    if (crereal)
		call adivkr (Memr[oRe], xnpts2, Memr[oRe], npts)
	    if (creimag)
		call adivkr (Memr[oIm], xnpts2, Memr[oIm], npts)
	}

	call ft_work_free (work, coord_shift, decenter)

	# Transpose the scratch file back into the output file.
	# If a file is not supposed to exist, delete it & reset the flag.

	if (crereal) {
	    call ft_transpose (Re_scr, FT_REPT(fto), len_blk)
	} else {
	    call ft_del_s (FT_REPT(fto))	# delete it
	    FT_REAL(fto) = NO
	}
	call ft_del_s (Re_scr)			# delete it

	if (creimag) {
	    call ft_transpose (Im_scr, FT_IMPT(fto), len_blk)
	} else {
	    call ft_del_s (FT_IMPT(fto))
	    FT_IMAG(fto) = NO
	}
	call ft_del_s (Im_scr)

	call mfree (work, TY_STRUCT)
end
