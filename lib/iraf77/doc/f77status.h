# IRAF77 -- Symbolic parameter definitions for F77/VOS interface


# Status return codes that match in SPP and Fortran includes

define	ER_EOF                 -2  # End-of-file indicator; same as IRAF
define	ER_OK                   0  # Successful return

# CL parameter access status messages : codes 1-9

define	ER_CLEOF               -1  # End of file on CL string
define	ER_LSTEOF              -1  # CL end of list
define	ER_CLNOTFND             1  # CL parameter not found
define	ER_CLBADTYP             2  # CL parameter has bad data type
define	ER_CLUNDF               3  # CL parameter is undefined
define	ER_CLPUT                5  # Error putting parameter into CL
define	ER_CLBADFELEM           6  # First element is bad (must be > 0)
define	ER_CLBADNELEMS          7  # Number of elements is bad (must be > 0)

# Image access status messages : codes 10-29

define	ER_IMOPOLD             10  # Error opening existing file
define	ER_IMOPNEW             11  # Error opening new file
define	ER_IMOPCOP             12  # Error opening new copy file
define	ER_IMBADNAXIS          13  # Invalid NAXIS parameter
define	ER_IMBADDIMEN          14  # Invalid axis length parameter
define	ER_IMBADTYPE           15  # Invalid image pixel data type
define	ER_IMOFFSET            16  # Error returning offset to image data
define	ER_IMBADACCMOD         17  # Invalid image access mode
define	ER_IMCLOS              18  # Error closing existing file
define	ER_IMNOTIMAG           19  # Not an image !
define	ER_IMBADSEC            20  # Bad section specification for image
define	ER_IMBADEXTN           21  # Bad extension for image
define	ER_IMILLSEC            22  # Illegal section on new image/new copy
define	ER_IMREAD              23  # Error reading pixel file
define	ER_IMWRITE             24  # Error writing to pixel file
define	ER_IMSTATUNKPAR        25  # Get unknown image parameter
define	ER_IMSETUNKPAR         26  # Set unknown image parameter

# Terminal I/O status messages : codes 30-39

define	ER_UPTOUT              30  # Error writing to STDOUT
define	ER_UPTLOG              31  # Error writing to the log file
define	ER_UPTBADDEST          32  # Invalid message destination
define	ER_UPTBADPRIO          33  # Invalid message priority

# Header I/O status messages : 40-59

define	ER_HDRPARNF            40  # Header parameter not found
define	ER_HDRPARTY            41  # Illegal data type for header parameter
define	ER_HDRNOSP             42  # Out of space in image header
define	ER_HDRPARRDF           43  # Illegal attempt to redefine parameter
define	ER_BADINPAT            45  # Bad input template (for parm tpt ??)
define	ER_HDBADFELEM          46  # Bad input first element
define	ER_HDBADNELEM          47  # Bad input number of elements
define	ER_HDNODEL             48  # Cannot delete standard keyword
define	ER_HDDELNXKW           49  # Cannot delete nonexistent keyword
define	ER_IMNEX               50  # Image does not exist
define	ER_IMDELETE            51  # Error deleting image
define	ER_IMRENAME            52  # Error renaming image
define	ER_HDRIMNF             53  # Illegal image header
define	ER_HDRIMDSCR           54  # Image descriptor mismatch

# Template processing status messages : codes 60-69

define	ER_IMTOPN              60  # Error opening template
define	ER_IMTGETNEXT          61  # Error getting next image from template
define	ER_IMTMAXPAT           62  # Too many patterns in template filename
define	ER_IMTILLPAT           63  # Illegal filename template
define	ER_IMTCLST             64  # Error closing template

# Graphics processing status messages : codes 70-99

define	ER_BADGRAPHDEV         70  # Bad graphics device
define	ER_BADGRAPHMODE        71  # Bad graphics mode
define	ER_GRAPHOPEN           72  # Error opening graphics device
define	ER_GRAPHCLOSE          73  # Error closing graphics device
define	ER_GRAPHCLEAR          74  # Error clearing graphics device
define	ER_GRAPHSETWND         75  # Error setting window
define	ER_GRAPHGETWND         76  # Error getting window
define	ER_GRAPHLINE           77  # Error plotting a polyline
define	ER_GRAPHMARKER         78  # Error plotting a sequence of markers
define	ER_GRAPHBADNPTS        79  # Invalid number of points
define	ER_GRAPHBADMARK        80  # Bad marker (size is incorrect)
define	ER_GRAPHBADAXES        81  # Bad axes specification
define	ER_GRAPHBADJUST        82  # Bad justification for test
define	ER_GRAPHBADSIZE        83  # Bad character size
define	ER_GRAPHBADORIENT      84  # Bad orientation (for text)
define	ER_GRAPHLABEL          85  # Error writing label
define	ER_GRAPHTEXT           86  # Error writing text
define	ER_GRAPHSETVP          87  # Error setting viewport
define	ER_GRAPHGETVP          88  # Error getting viewport
define	ER_GRAPHSETPAR         89  # Error setting GIO parameter
define	ER_GRAPHGETPAR         90  # Error getting GIO parameter
define	ER_GRAPHILLCODE        91  # Bad reset code (in reset function)
define	ER_GRAPHCURWR          92  # Error writing cursor
define	ER_GRAPHPLOT           93  # Error in high level plot routine
define	ER_GRAPHERRBAR         94  # Error in high level error bar routine
define	ER_GRAPHISTGM          95  # Error in high level histogram plot
define	ER_BADCURSOR           96  # Bad cursor call
define	ER_GRAPHSCALE          97  # Error in scale procedure

# Dynamic memory status messages : codes 100-109

define	ER_DMINVDTYP          100  # Invalid dynamic memory data type
define	ER_DYNMEMALC          101  # Error allocating dynamic memory
define	ER_FREEDYNMEM         102  # Error freeing dynamic memory

# Utility character/number conversion status messages : 110-119

define	ER_NOTNUMBER          110  # Not a number
define	ER_OVFNUMBER          111  # Overflow
define	ER_NOSPACE            112  # Field width not wide enough
