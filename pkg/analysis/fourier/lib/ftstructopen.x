include "../fourier.h"

# This file contains the following subroutines for allocating and
# deallocating memory for the FT structure for an input or output image:
#
# ft_struct_open
# ft_struct_close
#
# Phil Hodge, 15-Jan-1996  Subroutines created.

# ft_struct_open -- allocate memory for FT struct

procedure ft_struct_open (ft)

pointer ft		# o: FT structure
#--

begin
	# Allocate memory for the FT structure itself.
	call malloc (ft, LEN_FT, TY_STRUCT)

	# Allocate memory for the names of the real and imaginary parts.
	call malloc (FT_NAME_R_PTR(ft), SZ_LINE, TY_CHAR)
	call malloc (FT_NAME_I_PTR(ft), SZ_LINE, TY_CHAR)
end

# ft_struct_close -- free memory for FT struct

procedure ft_struct_close (ft)

pointer ft		# io: FT structure
#--

begin
	# Free memory for the names of the real and imaginary parts.
	call mfree (FT_NAME_I_PTR(ft), TY_CHAR)
	call mfree (FT_NAME_R_PTR(ft), TY_CHAR)

	# Free memory for the FT structure itself.
	call mfree (ft, TY_STRUCT)
end
