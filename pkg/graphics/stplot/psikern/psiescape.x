include <gset.h>
include <psiescape.h>
include "psi.h"

# Memory management
define	Sx		Memc[sx+$1-1]

#---------------------------------------------------------------------------
.help psi_escape 9Nov94 source
.ih
NAME
psi_escape -- Pass a device dependent instruction on to the kernel. 
.ih
USAGE
call psi_escape (fn, instruction, nwords)
.ih
ARGUMENTS
.ls fn (int)
The function code.
.le
.ls instruction (short[ARB])
The instrunction data words.
.le
.ls nwords (int)
The length of the instruction (number of entries in instruction).
.le
.ih
DESCRIPTION
The codes are defined in the include file "psiescape.h", see that file
for a description of the codes.  For the detailed workings of each code,
reference the code below and see the corresponding routine.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_escape (fn, instruction, nwords)

int     fn                      # function code
short   instruction[ARB]	# instruction data words
int     nwords                  # length of instruction

# Declarations.
int     i                       # Generic integer.

pointer sp                      # Stack pointer.
pointer	sx			# Generic string.
pointer tx                      # Text attribute pointer.

include "psi.com"

begin
        if (g_debug) {
            call eprintf ("psi_escape: Executing code %d with array length %d.\n")
            call pargi (fn)
            call pargi (nwords)
        }
        
	if (nwords <= 0)
	    return

	call smark (sp)
	
	# Set DCS indicator
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
	
        # Execute the appropriate routines.
        switch (fn) {
        case PS_CODE:
            # Convert the array to a string and output the string.
	    call salloc (sx, nwords, TY_CHAR)
	    call achtsc (instruction, Sx(1), nwords)
	    Sx(nwords+1) = EOS
	    call psk_out (Sx(1))

            
        case PS_IMAGE_RED_LUT, PS_IMAGE_GREEN_LUT, PS_IMAGE_BLUE_LUT,
						   PS_GR_RED_LUT, PS_GR_GREEN_LUT, PS_GR_BLUE_LUT:
            call psi_download_lut (fn, instruction, nwords)
            
        case PS_ROMAN_FONT, PS_GREEK_FONT, PS_ITALIC_FONT, PS_BOLD_FONT:
	    call salloc (sx, nwords, TY_CHAR)
            call achtsc (instruction, Sx(1), nwords)
	    Sx(nwords+1) = EOS
            switch (fn) {
            case PS_ROMAN_FONT:
                call psi_change_font (GT_ROMAN, Sx(1))
            case PS_GREEK_FONT:
                call psi_change_font (GT_GREEK, Sx(1))
            case PS_ITALIC_FONT:
                call psi_change_font (GT_ITALIC, Sx(1))
            case PS_BOLD_FONT:
                call psi_change_font (GT_BOLD, Sx(1))
            }
            
        case PS_VARIABLE_SPACE:
            # Get the text pointer and set the appropriate flag.
            tx = PSI_TXAP(g_kt)
            TX_VARIABLE(tx) = instruction[1]
            
        case PS_DASH:
            i = instruction[1]
            call psk_defi ("DASH", i)
            
        case PS_DOT:
            i = instruction[1]
            call psk_defi ("DOT", i)
            
        case PS_SPACE:
            i = instruction[1]
            call psk_defi ("SPACE", i)
            
        case PS_FILL_PATTERN:
            call smark(sp)
	    call salloc (sx, SZ_LINE, TY_CHAR)
            
            call sprintf (g_output, SZ_LINE, "/H%d <")
            call pargs (instruction[1])
            do i = 2, nwords {
                call sprintf (Sx(1), SZ_LINE, "%2.2x")
                call pargs (instruction[i])
                call strcat (Sx(1), g_output, SZ_LINE)
            }
            call strcat ("> /def", g_output, SZ_LINE)
            call psk_out (g_output)
            
            call sfree(sp)

	case PS_IMAGE_LUT:
	    call salloc (sx, nwords, TY_CHAR)
	    call achtsc (instruction, Sx(1), nwords)
	    Sx(nwords+1) = EOS
	    call psi_set_lut (PS_IMAGE, Sx(1))
            
	case PS_GRAPHICS_LUT:
	    call salloc (sx, nwords, TY_CHAR)
	    call achtsc (instruction, Sx(1), nwords)
	    Sx(nwords+1) = EOS
	    call psi_set_lut (PS_GRAPHIC, Sx(1))
            
        default:
	    if (g_debug) {
		call eprintf ("ERROR:psi_escape: Unknown escape sequence %d.\n")
		call pargi (fn)
	    }
        }

        call sfree (sp)
end
#---------------------------------------------------------------------------
# End of psi_escape
#---------------------------------------------------------------------------
