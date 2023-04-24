include <error.h>		# for EA_WARN
include	<imhdr.h>
include	"../fourier.h"
include	"../fterr.h"

# This file contains ft_open_files_i, ft_open_files_o, and ft_close_files.
#
# C. D. Biemesderfer, STScI, 9 Dec 1987  Original.
# Phil Hodge, 18-Apr-1990  ft_open_files_o:  set output data type to be real.
# Phil Hodge, 19-May-1993  Check dimension of input image in ft_open_files_i.
# Phil Hodge, 20-Jul-1993  Change file name convention; use getreal and
#			getimag to determine file names and to actually
#			specify which images are to be opened.
# Phil Hodge, 27-Jun-1994  Use iferr when calling load_ctstruct.
# Phil Hodge, 15-Jan-1996  Save names of real & imaginary parts in ft struct.

# FT_OPEN_FILES_I -- Open input files.
#
# The real and imaginary parts of the complex data are kept in separate images
# at this time.  The file names for these parts are related via a naming
# convention to the input and output file names specified by the user.  This
# convention is the obvious one, viz., the user-specified names with 'r' and
# 'i' appended.  (Note that this means that the user may specify file names
# that are only SZ_LINE-1 characters long.)  This routine also has to get
# world coordinate information from the input files and load them into the
# ft structure.
#
# As of July 1993, the convention is as follows.  If both getreal = true
# and getimag = true then "r" and "i" are appended to form the file names.
# If getreal = true but getimag = false or vice versa, then NOTHING is
# appended.  Also, if either flag is true the corresponding file MUST
# exist.  It used to be that if the flag was true but the file didn't
# exist, it was quietly ignored.

procedure ft_open_files_i (fti, infile, getreal, getimag)

pointer	fti		# i: FFT pointer; input image
char	infile[ARB]	# i: Input file name
bool	getreal		# i: Use input real part flag
bool	getimag		# i: Use input imaginary part flag
#--
pointer	immap(), irp, iip

begin
	# Apply naming convention.
	if (getreal && getimag) {
	    call ft_fname (infile, "r", FT_NAME_R(fti), SZ_LINE) # append "r"
	    call ft_fname (infile, "i", FT_NAME_I(fti), SZ_LINE) # append "i"
	} else if (getreal) {
	    call strcpy (infile, FT_NAME_R(fti), SZ_LINE)
	} else if (getimag) {
	    call strcpy (infile, FT_NAME_I(fti), SZ_LINE)
	} else {
	    call error (1, "no input:  getreal=no and getimag=no")
	}

	# Open input real data.
	if (getreal) {
	    irp = immap (FT_NAME_R(fti), READ_ONLY, 0)
	    FT_REAL(fti) = YES
	    FT_REPT(fti) = irp

	    if (IM_NDIM(irp) > FT_MAX_NDIM) {
		call imunmap (irp)
		call ft_error (FT_BADNAXIS, FT_FATAL)
	    }
	} else {
	    call strcpy ("", FT_NAME_R(fti), SZ_LINE)
	    FT_REAL(fti) = NO
	    FT_REPT(fti) = NULL
	}

	# Open input imaginary data.
	if (getimag) {
	    iip = immap (FT_NAME_I(fti), READ_ONLY, 0)
	    FT_IMAG(fti) = YES
	    FT_IMPT(fti) = iip

	    if (IM_NDIM(iip) > FT_MAX_NDIM) {
		call imunmap (iip)
		call ft_error (FT_BADNAXIS, FT_FATAL)
	    }
	} else {
	    call strcpy ("", FT_NAME_I(fti), SZ_LINE)
	    FT_IMAG(fti) = NO
	    FT_IMPT(fti) = NULL
	}

	# Preferentially use input real file as template.  Select imaginary
	# file only if there is no real part specified.  If neither one is
	# specified, there is no input, and that's an error.

	if (FT_REAL(fti) == YES)
	    FT_IMAGE(fti) = irp
	else if (FT_IMAG(fti) == YES)
	    FT_IMAGE(fti) = iip
	else
	    call ft_error (FT_NOINPUT, FT_FATAL)

	# Get size and coordinate information.
	iferr {
	    call load_ctstruct (fti, FT_IMAGE(fti))
	} then {
	    call erract (EA_WARN)
	    call ft_wcs_0 (fti)		# assign default WCS
	}
end

# FT_OPEN_FILES_O -- Open output files.
#
# The real and imaginary parts of the complex data are kept in separate images
# at this time.  The file names for these parts are related via a naming
# convention to the input and output file names specified by the user.  This
# convention is the obvious one, viz., the user-specified names with 'r' and
# 'i' appended.  (Note that this means that the user may specify file names
# that are only SZ_LINE-1 characters long.)

procedure ft_open_files_o (fti, fto, outfile, crereal, creimag)

pointer fti		# i: FFT pointer; input image
pointer fto		# i: FFT pointer; output image
char	outfile[ARB]	# i: Output file name
bool	crereal		# i: Use output real part flag
bool	creimag		# i: Use output imaginary part flag
#--
pointer	immap(), orp, oip

begin
	# Apply naming convention.
	if (crereal && creimag) {
	    call ft_fname (outfile, "r", FT_NAME_R(fto), SZ_LINE) # append "r"
	    call ft_fname (outfile, "i", FT_NAME_I(fto), SZ_LINE) # append "i"
	} else if (crereal) {
	    call strcpy (outfile, FT_NAME_R(fto), SZ_LINE)
	} else if (creimag) {
	    call strcpy (outfile, FT_NAME_I(fto), SZ_LINE)
	} else {
	    call error (1, "no output:  crereal=no and creimag=no")
	}

	# Open output real data.
	if (crereal) {
	    orp = immap (FT_NAME_R(fto), NEW_COPY, FT_IMAGE(fti))
	    FT_REAL(fto) = YES
	    FT_REPT(fto) = orp
	    IM_PIXTYPE(orp) = TY_REAL		# in case input is int or short
	} else {
	    call strcpy ("", FT_NAME_R(fto), SZ_LINE)
	    FT_REAL(fto) = NO
	    FT_REPT(fto) = NULL
	}

	# Open output imaginary data.
	if (creimag) {
	    oip = immap (FT_NAME_I(fto), NEW_COPY, FT_IMAGE(fti))
	    FT_IMAG(fto) = YES
	    FT_IMPT(fto) = oip
	    IM_PIXTYPE(oip) = TY_REAL
	} else {
	    call strcpy ("", FT_NAME_I(fto), SZ_LINE)
	    FT_IMAG(fto) = NO
	    FT_IMPT(fto) = NULL
	}

	if (FT_REAL(fto) == YES)
	    FT_IMAGE(fto) = orp
	else if (FT_IMAG(fto) == YES)
	    FT_IMAGE(fto) = oip
	else
	    call ft_error (FT_NOOUTPUT, FT_FATAL)

	FT_NAXIS(fto) = FT_NAXIS(fti)
end

# FT_CLOSE_FILES -- Close real & imaginary parts

procedure ft_close_files (ft)

pointer ft			# i: FFT pointer
#--

begin
	# Unmap images
	if (FT_IMAG(ft) == YES)
	    call imunmap (FT_IMPT(ft))

	if (FT_REAL(ft) == YES)
	    call imunmap (FT_REPT(ft))
end
