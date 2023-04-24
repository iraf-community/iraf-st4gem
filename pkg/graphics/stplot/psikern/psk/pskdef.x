include "psk.h"

#---------------------------------------------------------------------------
.help psk_def[birsl] 1May92 plot
.ih
NAME
psk_def[birsl] -- Setup PostScript definitions.
.ih
USAGE
call psk_def[birsl] (variable, value)
.ih
ARGUMENTS
.ls variable (char[ARB])
The name of the PostScript variable to define.
.le
.ls value (bool|int|real|char[ARB]|char[ARB])
The initial value to give the variable.
.le
.ih
SEE ALSO
psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_defb (variable, value)

char variable[ARB]  # I:  Name of the PostScrip variable to define.
bool value          # I:  The value to give it.

# Declarations.
include "psk.com"

begin
        
        if (ps_debug) {
            call eprintf ("Defining bool variable %s with value %b.\n")
            call pargstr (variable)
            call pargb (value)
        }
        
        if (value) {
            call sprintf (ps_output, MAX_CHAR, "/%s true def")
            call pargstr (variable)
        } else {
            call sprintf (ps_output, MAX_CHAR, "/%s false def")
            call pargstr (variable)
        }
        call psk_out (ps_output)
        
end
#---------------------------------------------------------------------------
# End of psk_defb
#---------------------------------------------------------------------------
procedure psk_defi (variable, value)

char variable[ARB]  # I:  Name of the PostScrip variable to define.
int  value          # I:  The value to give it.

# Declarations.
include "psk.com"

begin
        
        if (ps_debug) {
            call eprintf ("Defining int variable %s with value %d.\n")
            call pargstr (variable)
            call pargi (value)
        }
        
        call sprintf (ps_output, MAX_CHAR, "/%s %d def")
        call pargstr (variable)
        call pargi (value)
        call psk_out (ps_output)
        
end
#---------------------------------------------------------------------------
# End of psk_defi
#---------------------------------------------------------------------------
procedure psk_defr (variable, value)

char variable[ARB]  # I:  Name of the PostScrip variable to define.
real value          # I:  The value to give it.

# Declarations.
include "psk.com"

begin
        
        if (ps_debug) {
            call eprintf ("Defining real variable %s with value %g.\n")
            call pargstr (variable)
            call pargr (value)
        }
        
        call sprintf (ps_output, MAX_CHAR, "/%s %g def")
        call pargstr (variable)
        call pargr (value)
        call psk_out (ps_output)
        
end
#---------------------------------------------------------------------------
# End of psk_defr
#---------------------------------------------------------------------------
procedure psk_defs (variable, value)

char variable[ARB]  # I:  Name of the PostScrip variable to define.
char value[ARB]     # I:  The value to give it.

# Declarations.
pointer sp   # Stack pointer
pointer tmp  # Temporary string.

include "psk.com"

# Function prototypes
int strlen()

begin
        
        call smark (sp)
        call salloc (tmp, 3 * strlen (value), TY_CHAR)
        
        if (ps_debug) {
            call eprintf ("Defining string variable %s with value %s.\n")
            call pargstr (variable)
            call pargstr (value)
        }
        
        call psi_mkpsstring (value, Memc[tmp], 3 * strlen (value))
        call sprintf (ps_output, MAX_CHAR, "/%s %s def")
        call pargstr (variable)
        call pargstr (Memc[tmp])
        call psk_out (ps_output)
        
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of psk_defs
#---------------------------------------------------------------------------
procedure psk_defl (variable, value)

char variable[ARB]  # I:  Name of the PostScrip variable to define.
char value[ARB]     # I:  The value to give it.

include "psk.com"

begin
        
        if (ps_debug) {
            call eprintf ("Defining literal variable %s with value %s.\n")
            call pargstr (variable)
            call pargstr (value)
        }
        
        call sprintf (ps_output, MAX_CHAR, "/%s /%s def")
        call pargstr (variable)
        call pargstr (value)
        call psk_out (ps_output)
        
end
#---------------------------------------------------------------------------
# End of psk_defl
#---------------------------------------------------------------------------
