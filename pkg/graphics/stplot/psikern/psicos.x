# Define the length of the output string.
define LENGTH 5

#---------------------------------------------------------------------------
.help psi_compress_coords 1May92 plot
.ih
NAME
psi_compress_coords -- Produce a string representing compressed coordinates.
.ih
USAGE
<stringpointer> = psi_compress_coords (x, y)
.ih
ARGUMENTS
.ls x, y (int)
The coordinates to be represented in a compressed string. See description
below.
.le
.ih
RETURNS
.ls <stringpointer>
A pointer to a string representing the compressed coordinates.
.le
.ih
DESCRIPTION
In order to reduce the number of bytes sent to a PostScript device, this
routine produces a string that represents some compressed format for
the two coordinates.  The format is defined by the PostScript procedure
defined in psikern_prolog.ps (or whatever prolog file is being used).
NOTE: The coordinates must already have been properly scaled for whatever
compression routine is used.  At the moment, both coordinates must satisfy
the following relation: 0 <= x,y <= 4095.

When the calling routine is done with the string, it may use the mfree
call to deallocate the memory.

NOTE: This routine comes from the IRAF sgi2uapl SGI translator.
.ih
SEE ALSO
t_psikern, psikern_prolog, psi_draw, psi_move
.endhelp
#---------------------------------------------------------------------------

pointer procedure psi_compress_coords (x, y)

int x, y  # I:  The coordinates to be written to a string.

# Declarations.
int   n         # The encoded value.

pointer output  # The pointer to the output string.

# Function prototypes
int andi(), ori(), shifti()

begin

        # Allocate the string memory.
        call malloc (output, LENGTH, TY_CHAR)
        
        # Convert the X coordinate.
        n = andi (shifti (x, -6), 77b)
        if (n != 77b)
            n = ori (n, 100b)
        Memc[output] = n
        n = andi (x, 77b)
        if (n != 77b)
            n = ori (n, 100b)
        Memc[output+1] = n
        
        # Convert the Y coordinate.
        n = andi (shifti (y, -6), 77b)
        if (n != 77b)
            n = ori (n, 100b)
        Memc[output+2] = n
        n = andi (y, 77b)
        if (n != 77b)
            n = ori (n, 100b)
        Memc[output+3] = n
        
        Memc[output+4] = EOS
        
        # Return the string.
        return (output)
        
end
#---------------------------------------------------------------------------
# End of psi_compress_coords
#---------------------------------------------------------------------------
