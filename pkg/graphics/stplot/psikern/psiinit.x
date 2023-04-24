include <ctype.h>
include <gki.h>
include <mach.h>
include "psi.h"

# Define the defaults for dashs and dots
define DEFAULT_DASH_SIZE   400
define DEFAULT_DOT_SIZE    40
define DEFAULT_SPACE_SIZE  200

# Defaults for line widths.
define DEFAULT_LWSLOPE 10.

#---------------------------------------------------------------------------
.help psi_init 4Nov94 source
.ih
NAME
psi_init -- Initialize the gkt data structures.
.ih
USAGE
call psi_init (tty, devname)
.ih
ARGUMENTS
.ls tty (pointer)
Graphcap file descriptor
.le
.ls devname (char[ARB])
Device name.
.le
.ih
DESCRIPTION
Initialize the gkt data structures from the graphcap entry
for the device.  Called once, at OPENWS time, with the TTY pointer already
set in the common.  The companion routine PSI_RESET initializes the attribute
packets when the frame is flushed.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_init (tty, devname)

pointer tty                     # graphcap descriptor
char    devname[ARB]            # device name

# Declarations.
real char_height     # Character sizes
real char_size
real char_width

int i                # Generic.
int maxch            # Maximum size of string.

pointer nextch       # Pointer to next character.
pointer prolog_file  # The prolog file name.
pointer sp           # Stack pointer.
pointer	xp		# Generic.

include "psi.com"
include "psiparams.com"

# Function prototypes.
real    ttygetr()
int     btoi(), errcode(), gstrcpy(), open(), strlen()
int	tty_find_capability()
int     ttygeti(), ttygets()
bool    ttygetb()

begin
        
        # Allocate the gkt descriptor and the string buffer.
        if (g_kt == NULL) {
            call calloc (g_kt, LEN_PSI, TY_STRUCT)
            call malloc (PSI_SBUF(g_kt), SZ_SBUF, TY_CHAR)
        }
        
        # Allocate memory.
        call smark (sp)
        call salloc (prolog_file, SZ_FNAME, TY_CHAR)
        
        # Init string buffer parameters.  The first char of the string buffer
        # is reserved as a null string, used for graphcap control strings
        # omitted from the graphcap entry for the device.
        PSI_SZSBUF(g_kt) = SZ_SBUF
        PSI_NEXTCH(g_kt) = PSI_SBUF(g_kt) + 1
        Memc[PSI_SBUF(g_kt)] = EOS
        
        # See if we should be debugging.
        g_debug = ttygetb (tty, "DB")
        
        # Determine device scale and offset from graphcap entries.
        PSI_PORTRAIT(g_kt) = ttygetb (tty, "PT")
        
        PSI_XOFF(g_kt) = ttygetr (tty, "XO")
        PSI_YOFF(g_kt) = ttygetr (tty, "YO")
        
        PSI_XSIZE(g_kt) = ttygetr (tty, "xs")
        PSI_YSIZE(g_kt) = ttygetr (tty, "ys")
        
        # Get the line width and increment.
        PSI_LWORIGIN(g_kt) = GKI_MAXNDC * ttygetr (tty, "PW")
        PSI_LWSLOPE(g_kt) = ttygetr (tty, "PI")
        
        # Initialize the scale factor.
        PSI_GKI2OUT(g_kt) = real (PS_OUT_RESOLUTION) / REAL (GKI_MAXNDC)
        
        # Initialize the character scaling parameters, required for text
        # generation.  The heights are given in NDC units in the graphcap
        # file, which we convert to GKI units.  Estimated values are
        # supplied if the parameters are missing in the graphcap entry.
        char_height = ttygetr (tty, "ch")
        char_height = char_height * GKI_MAXNDC
        char_width = ttygetr (tty, "cw")
        char_width = char_width * GKI_MAXNDC
        
        # If the device has a set of discreet character sizes, get the
        # size of each by fetching the parameter "tN", where the N is
        # a digit specifying the text size index.  Compute the height and
        # width of each size character from the "ch" and "cw" parameters
        # and the relative scale of character size I.
        # NOTE: This is, in general, not applicable to PostScript devices,
        # but there is no reason not to allow this capability.
        PSI_NCHARSIZES(g_kt) = min (MAX_CHARSIZES, ttygeti (tty, "th"))
        nextch = PSI_NEXTCH(g_kt)
        
        if (PSI_NCHARSIZES(g_kt) <= 0) {
            PSI_NCHARSIZES(g_kt) = 1
            PSI_CHARSIZE(g_kt,1) = 1.0
            PSI_CHARHEIGHT(g_kt,1) = char_height
            PSI_CHARWIDTH(g_kt,1)  = char_width
        } else {
            Memc[nextch+2] = EOS
            for (i=1;  i <= PSI_NCHARSIZES(g_kt);  i=i+1) {
                Memc[nextch] = 't'
                Memc[nextch+1] = TO_DIGIT(i)
                char_size = ttygetr (tty, Memc[nextch])
                PSI_CHARSIZE(g_kt,i)   = char_size
                PSI_CHARHEIGHT(g_kt,i) = char_height * char_size
                PSI_CHARWIDTH(g_kt,i)  = char_width  * char_size
            }
        }
        
        # Get the Dash, Dot, and Space sizes, if defined.
        PSI_DASH(g_kt) = GKI_MAXNDC * ttygetr (tty, "TD")
        if (PSI_DASH(g_kt) == 0)
            PSI_DASH(g_kt) = DEFAULT_DASH_SIZE
        PSI_DASH(g_kt) = PSI_DASH(g_kt) * PSI_GKI2OUT(g_kt)
        
        PSI_DOT(g_kt) = GKI_MAXNDC * ttygetr (tty, "TP")
        if (PSI_DOT(g_kt) == 0)
            PSI_DOT(g_kt) = DEFAULT_DOT_SIZE
        PSI_DOT(g_kt) = PSI_DOT(g_kt) * PSI_GKI2OUT(g_kt)
        
        PSI_SPACE(g_kt) = GKI_MAXNDC * ttygetr (tty, "TS")
        if (PSI_SPACE(g_kt) == 0)
            PSI_SPACE(g_kt) = DEFAULT_SPACE_SIZE
        PSI_SPACE(g_kt) = PSI_SPACE(g_kt) * PSI_GKI2OUT(g_kt)
        
        # Get the prolog file name and open the file.
        if (ttygets (tty, "IF", Memc[prolog_file], SZ_FNAME) <= 0)
            call error (1, "psi_init: Could not find prolog file graphcap entry IF!")
        iferr (PSI_PROLOG(g_kt) = open (Memc[prolog_file], READ_ONLY, TEXT_FILE)) {
            i = errcode()
            call eprintf ("psi_init: Could not open prolog file %s.\n")
            call pargstr (Memc[prolog_file])
            call erract (i)
        }
        
        # Get the name of the fonts, if any are specified.  
        if (strlen (roman_font) == 0)
            if (ttygets (tty, "FR", roman_font, SZ_LINE) <= 0)
                call strcpy ("Times-Roman", roman_font, SZ_LINE)
        
        if (strlen (greek_font) == 0)
            if (ttygets (tty, "FG", greek_font, SZ_LINE) <= 0)
                call strcpy ("Symbol", greek_font, SZ_LINE)
        
        if (strlen (bold_font) == 0)
            if (ttygets (tty, "FB", bold_font, SZ_LINE) <= 0)
                call strcpy ("Times-Bold", bold_font, SZ_LINE)
        
        if (strlen (italic_font) == 0)
            if (ttygets (tty, "FI", italic_font, SZ_LINE) <= 0)
                call strcpy ("Times-Italic", italic_font, SZ_LINE)

	# Get LUT names.
        if (strlen (lut_gr) == 0)
            i = ttygets (tty, "LG", lut_gr, SZ_LINE)
        if (strlen (lut_im) == 0)
            i = ttygets (tty, "LI", lut_im, SZ_LINE)
	
        # Get the device resolution.
        PSI_XRES(g_kt)       =       ttygeti (tty, "xr")
        PSI_YRES(g_kt)       =       ttygeti (tty, "yr")
        PSI_ZRES(g_kt)       =       ttygeti (tty, "zr")
        
        # Determine whether the default text is variable or mono-spaced.
        PSI_DEFVAR(g_kt) =  btoi (!ttygetb (tty, "MO"))
        if (!IS_INDEFI(proportional))
            PSI_DEFVAR(g_kt) = proportional
        
        # Get the number of frames per metafile.
        g_maxframes = ttygeti (tty, "MF")
        if (strlen (out_file) > 0)
            g_maxframes = MAX_INT
        
        # Initialize the input parameters.
        PSI_CURSOR(g_kt)     = 1

	# Determine whether fill areas and cell arrays will be supported.
	PSI_DOCELL(g_kt) = btoi (ttygetb (tty, "ca"))
	PSI_DOFILL(g_kt) = btoi (ttygetb (tty, "fa"))

        # Save the device string in the descriptor.
        nextch = PSI_NEXTCH(g_kt)
        PSI_DEVNAME(g_kt) = nextch
        maxch = PSI_SBUF(g_kt) + SZ_SBUF - nextch + 1
        nextch = nextch + gstrcpy (devname, Memc[nextch], maxch) + 1
        
        PSI_NEXTCH(g_kt) = nextch

	# Determine whether the background should be painted, and
	# what color it should be.
	PSI_BACK(g_kt) = INDEFI
	if (tty_find_capability (tty, "PB", xp) == YES)
	    PSI_BACK(g_kt) = ttygeti (tty, "PB")

        # Release memory.
        call sfree (sp)
end
#---------------------------------------------------------------------------
# End of psi_init
#---------------------------------------------------------------------------
