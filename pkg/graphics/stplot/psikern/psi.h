.help psi.h 1May92 plot
.ih
NAME
psi global definitions.
.ih
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

define  MAX_CHARSIZES   10                      # max discreet device char sizes
define  SZ_SBUF         1024                    # initial string buffer size
define  SZ_GDEVICE      31                      # maxsize forced device name

# Define the output resolution.
define PS_OUT_RESOLUTION 4095

# Define names for the color lookup tables.
define PS_IMAGE    1
define PS_GRAPHIC  2

# How many line segments can be drawn without saving them to the page.
define  PS_MAX_SEGMENTS         500

# Define names for using USER specified HARD fonts
#	PS_SETUP signals to output the '/PF /<font-name> def' string to 
#		output postscript
#	PS_USER corresponds to the escape sequence used in the *label 
#		commands in IGI to use the /PF font in the Postscript output
# These definitions extend the current 'GT_ROMAN,...' definitions in <gset.h>
define	PS_USER		100
define	PS_SETUP	101

# Define names for the seperate colors of lookup tables.
define PS_RED    1
define PS_GREEN  2
define PS_BLUE   3

# What should be considered the default color.
define  DEF_COLOR       1
define  DEF_MAX_COLOR   15

# The PSI state/device descriptor.

define  LEN_PSI         101

define  PSI_SBUF        Memi[$1]                # string buffer
define  PSI_SZSBUF      Memi[$1+1]              # size of string buffer
define  PSI_NEXTCH      Memi[$1+2]              # next char pos in string buf
define  PSI_NCHARSIZES  Memi[$1+3]              # number of character sizes
define  PSI_DEFVAR      Memi[$1+4]              # Default text spacing.
define  PSI_FONTLIST    Memi[$1+5]              # Pointer to the font list.
define  PSI_XRES        Memi[$1+6]              # device resolution in X
define  PSI_YRES        Memi[$1+7]              # device resolution in Y
define  PSI_ZRES        Memi[$1+8]              # device resolution in Z
define  PSI_FILLSTYLE   Memi[$1+9]              # number of fill styles
define  PSI_LWSLOPE     Memr[P2R($1+10)]        # Difference in line widths in
                                                # GKI units.
define  PSI_LWORIGIN    Memr[P2R($1+11)]        # Line widht in GKI units of
                                                # a line of width 1.
define  PSI_CURSOR      Memi[$1+12]             # last cursor accessed
define  PSI_COLOR       Memi[$1+13]             # last color set
define  PSI_TXSIZE      Memi[$1+14]             # last text size set
define  PSI_TXFONT      Memi[$1+15]             # last text font set
define  PSI_TYPE        Memi[$1+16]             # last line type set
define  PSI_WIDTH       Memi[$1+17]             # last line width set
define  PSI_DEVNAME     Memi[$1+18]             # name of open device
define  PSI_PORTRAIT    Memb[$1+19]             # True to print out in portait.
define  PSI_XSIZE       Memr[P2R($1+20)]        # Page size along X (meters)
define  PSI_YSIZE       Memr[P2R($1+21)]        # Page size along Y (meters)
define  PSI_XOFF        Memr[P2R($1+22)]        # Page offset along X (meters).
define  PSI_YOFF        Memr[P2R($1+23)]        # Page offset along Y (meters).
define  PSI_PROLOG      Memi[$1+24]             # PostScrip prolog file desc.
define  PSI_UP          Memi[$1+25]             # Last text angle.
define  PSI_GKI2OUT     Memr[P2R($1+26)]        # Scale conversion from GKI
                                                # to a reduced output space.
define  PSI_DASH        Memi[$1+27]             # Size of a dash in output space
define  PSI_DOT         Memi[$1+28]             # Size of dot in output space.
define  PSI_SPACE       Memi[$1+29]             # Size of space in output space.
define  PSI_PATH        Memi[$1+30]             # Last text path set.
define  PSI_HJUSTIFY    Memi[$1+31]             # Last horizontal justification.
define  PSI_VJUSTIFY    Memi[$1+32]             # Last verticle justification.
define  PSI_VARIABLE    Memi[$1+33]             # Last state of variable-spaced.
define  PSI_SPACING     Memr[P2R($1+34)]        # Space relative to height.
define	PSI_BACK	Memi[$1+35]		# Background color.
define	PSI_DOCELL	Memi[$1+36]		# Draw cell arrays?
define	PSI_DOFILL	Memi[$1+37]		# Do area fills?
        # extra space
define  PSI_CHARHEIGHT  Memi[$1+50+$2-1]        # character height
define  PSI_CHARWIDTH   Memi[$1+60+$2-1]        # character width
define  PSI_CHARSIZE    Memr[P2R($1+70+$2-1)]   # text sizes permitted
define  PSI_PLAP        ($1+80)                 # polyline attributes
define  PSI_PMAP        ($1+84)                 # polymarker attributes
define  PSI_FAAP        ($1+88)                 # fill area attributes
define  PSI_TXAP        ($1+91)                 # default text attributes

# Substructure definitions.

define  LEN_PL          4
define  PL_STATE        Memi[$1]                # polyline attributes
define  PL_LTYPE        Memi[$1+1]
define  PL_WIDTH        Memi[$1+2]
define  PL_COLOR        Memi[$1+3]

define  LEN_PM          4
define  PM_STATE        Memi[$1]                # polymarker attributes
define  PM_LTYPE        Memi[$1+1]
define  PM_WIDTH        Memi[$1+2]
define  PM_COLOR        Memi[$1+3]

define  LEN_FA          3                       # fill area attributes
define  FA_STATE        Memi[$1]
define  FA_STYLE        Memi[$1+1]
define  FA_COLOR        Memi[$1+2]

define  LEN_TX          11                      # text attributes
define  TX_STATE        Memi[$1]
define  TX_UP           Memi[$1+1]
define  TX_SIZE         Memi[$1+2]
define  TX_PATH         Memi[$1+3]
define  TX_SPACING      Memr[P2R($1+4)]
define  TX_HJUSTIFY     Memi[$1+5]
define  TX_VJUSTIFY     Memi[$1+6]
define  TX_FONT         Memi[$1+7]
define  TX_QUALITY      Memi[$1+8]
define  TX_COLOR        Memi[$1+9]
define  TX_VARIABLE     Memi[$1+10]

#---------------------------------------------------------------------------
# End of psi.h
#---------------------------------------------------------------------------
