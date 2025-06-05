C IRAF77 -- Fortran-77 include file for the F77/VOS interface
C Modified: D.L.Ball, 6-OCT-87
C             Remove tabs, change to uppercase
C           D.L.Ball, 15-DEC-1987
C             Changed '777O' to 511 for GRSETA
C             Comment out TYTEXT from datatypes; not used
C           C.D.Biemesderfer, 26-JAN-1988
C             Remove MEMC
C             Add US* status return parameters
C------------------------------------------------------------------------------
C Get IRAF MEM common into main program.
C
      LOGICAL          MEMB(1)
      INTEGER*2        MEMS(1)
      INTEGER*4        MEMI(1)
      INTEGER*4        MEML(1)
      REAL             MEMR(1)
      DOUBLE PRECISION MEMD(1)
      COMPLEX          MEMX(1)
      EQUIVALENCE (MEMB, MEMS, MEMI, MEML, MEMR, MEMD, MEMX)
      COMMON /MEM/ MEMD
C------------------------------------------------------------------------------
C Include F77/VOS parameter equivalents for stuff in hlib$iraf.h.
C
C                                       End-of-file code
      INTEGER    EOF
      PARAMETER (EOF = -2)
C                                       IRAF definition for NO
      INTEGER    NO
      PARAMETER (NO = 0)
C                                       IRAF definition for YES
      INTEGER    YES
      PARAMETER (YES = 1)
C                                       File access modes, etc.:
C                                       ------------------------
C                                       Read only
      INTEGER    RDONLY
      PARAMETER (RDONLY = 1)
C                                       Read/Write
      INTEGER    RDWRIT
      PARAMETER (RDWRIT = 2)
C                                       Write only
      INTEGER    WRONLY
      PARAMETER (WRONLY = 3)
C                                       Append
      INTEGER    APPEND
      PARAMETER (APPEND = 4)
C                                       New file
      INTEGER    NEWFIL
      PARAMETER (NEWFIL = 5)
C                                       Temporary file
      INTEGER    TMPFIL
      PARAMETER (TMPFIL = 6)
C                                       New copy image
      INTEGER    NEWCPY
      PARAMETER (NEWCPY = 7)
C                                       New image
      INTEGER    NEWIMG
      PARAMETER (NEWIMG = 5)
C                                       Data types:
C                                       -----------
C                                       Boolean
      INTEGER    TYBOOL
      PARAMETER (TYBOOL = 1)
C                                       Short integer
      INTEGER    TYSHOR
      PARAMETER (TYSHOR = 3)
C                                       Integer
      INTEGER    TYINT
      PARAMETER (TYINT = 4)
C                                       Long integer
      INTEGER    TYLONG
      PARAMETER (TYLONG = 5)
C                                       Single-precision real
      INTEGER    TYREAL
      PARAMETER (TYREAL = 6)
C                                       Double-precision real
      INTEGER    TYDOUB
      PARAMETER (TYDOUB = 7)
C                                       Complex
      INTEGER    TYCPLX
      PARAMETER (TYCPLX = 8)
C                                       Unsigned short integer
      INTEGER    TYUSHT
      PARAMETER (TYUSHT = 11)
C                                       Unsigned byte
      INTEGER    TYUBYT
      PARAMETER (TYUBYT = 12)
C
C TYTEXT is a special code for the IRAF77 interface; it is not in the VOS.
C
C The following is NOT the value to use when defining a character-type column
C for a table; use -N, where N is the maximum number of characters in an
C element of that column.
C
C NOT PRESENTLY USED IN THE INTERFACE; commented out 12/15/87, DB
C
C                                       Text string
C      INTEGER    TYTEXT
C      PARAMETER (TYTEXT = 13)
C                                       Indefinites:
C                                       ------------
C                                       Undefined short integer
      INTEGER    INDEFS
      PARAMETER (INDEFS = -32768)
C                                       Undefined long integer
      INTEGER    INDEFL
      PARAMETER (INDEFL = -2147483647)
C                                       Undefined integer
      INTEGER    INDEFI
      PARAMETER (INDEFI = INDEFL)
C                                       Undefined single-precision real
      REAL       INDEFR
      PARAMETER (INDEFR = 1.6E38)
C                                       Undefined double-precision real
      DOUBLE PRECISION INDEFD
      PARAMETER (INDEFD = 1.6D38)
C                                       Undefined complex
      COMPLEX    INDEFX
      PARAMETER (INDEFX = (INDEFR, INDEFR))
C                                       Undefined generic type
      REAL       INDEF
      PARAMETER (INDEF = INDEFR)
C------------------------------------------------------------------------------
C Here are F77/VOS parameter equivalents for stuff in lib$imset.h.
C
C Codes for image processing parameters accessible by F77/VOS users.  Be
C careful; you could get yourself in a heap of trouble with some of these
C if you don't know what you are doing !
C
C (NOTE:  Not all of the defined values in lib$imset.h are included here; if
C programmers desire more flexibility, then coding in SPP is more appropriate.)
C
C					maximum number of dimensions
      INTEGER	 MAXDIM
      PARAMETER	(MAXDIM = 7)
C                                       Width of boundary region in pixels
      INTEGER    IMNBND
      PARAMETER (IMNBND = 4)
C                                       Types of boundary extension (see below)
      INTEGER    IMTYBN
      PARAMETER (IMTYBN = 5)
C                                       Set bad pixel to undefined ?
      INTEGER    IMFBPX
      PARAMETER (IMFBPX = 6)
C                                       Update header at close ?
      INTEGER    IMUPHD
      PARAMETER (IMUPHD = 8)
C                                       Set const value for out-of-bounds pixels
      INTEGER    IMBPXV
      PARAMETER (IMBPXV = 9)
C                                       Name of open image section
      INTEGER    IMIMNM
      PARAMETER (IMIMNM = 13)
C                                       Boundary extension types:
C                                       -------------------------
C                                       Return constant if out of bounds
      INTEGER    BTCONS
      PARAMETER (BTCONS = 1)
C                                       Return nearest boundary pixel
      INTEGER    BTNEAR
      PARAMETER (BTNEAR = 2)
C                                       Reflect back into image
      INTEGER    BTREFL
      PARAMETER (BTREFL = 3)
C                                       Wrap around to other side
      INTEGER    BTWRAP
      PARAMETER (BTWRAP = 4)
C                                       Project about boundary
      INTEGER    BTPROJ
      PARAMETER (BTPROJ = 5)
C------------------------------------------------------------------------------
C Parameter definitions for F77/VOS interface arguments.
C
C                                       UMSPUT destinations:
C                                       --------------------
C                                       Write message on stdout
      INTEGER    STDOUT
      PARAMETER (STDOUT = 1)
C                                       Write message on stderr
      INTEGER    STDERR
      PARAMETER (STDERR = 2)
C                                       Write message on user logfile
      INTEGER    USRLOG
      PARAMETER (USRLOG = 4)
C                                       Write message on system logfile
      INTEGER    SYSLOG
      PARAMETER (SYSLOG = 8)
C
C                                       UHDAS header parm types:
C                                       ------------------------
C                                       Global header parameter
      INTEGER    GENHDR
      PARAMETER (GENHDR = 0)
C                                       Image-specific header parameter
      INTEGER    IMSPEC
      PARAMETER (IMSPEC = 1)
C------------------------------------------------------------------------------
C Definitions for parameters relevant to table I/O.
C
C                                       These may be set by UTPPTI and/or
C                                         read by UTPGTI:
C                                       ---------------------------------
C                                       Length of row (unit: size of real)
      INTEGER    TBRLEN
      PARAMETER (TBRLEN = 1)
C                                       Increase row length
      INTEGER    TBIRLN
      PARAMETER (TBIRLN = 2)
C                                       Number of rows to allocate
      INTEGER    TBALLR
      PARAMETER (TBALLR = 3)
C                                       Increase alloc num of rows
      INTEGER    TBIALR
      PARAMETER (TBIALR = 4)
C                                       Which type of table (row or column) ?
      INTEGER    TBWTYP
      PARAMETER (TBWTYP = 5)
C                                       Maximum number of user parameters
      INTEGER    TBMXPR
      PARAMETER (TBMXPR = 6)
C                                       Maximum number of columns
      INTEGER    TBMXCL
      PARAMETER (TBMXCL = 7)
C                                       Type: row-ordered table
      INTEGER    TBTYPR
      PARAMETER (TBTYPR = 11)
C                                       Type: column-ordered table
      INTEGER    TBTYPC
      PARAMETER (TBTYPC = 12)
C                                       These may be read by UTPGTI but
C                                         may not be set:
C                                       -------------------------------
C                                       Number of rows written to
      INTEGER    TBNROW
      PARAMETER (TBNROW = 21)
C                                       Number of columns defined
      INTEGER    TBNCOL
      PARAMETER (TBNCOL = 22)
C                                       Amount of row used (unit: size of real)
      INTEGER    TBRUSD
      PARAMETER (TBRUSD = 23)
C                                       Number of user parameters
      INTEGER    TBNPAR
      PARAMETER (TBNPAR = 24)
C------------------------------------------------------------------------------
C Definitions for parameters relevant to graphics.
C
C                                       General graphics parameters:
C                                       ----------------------------
C                                       WCS number
      INTEGER    GRPWCS
      PARAMETER (GRPWCS = 3)
C                                       X WCS transformation
      INTEGER    GXTRAN
      PARAMETER (GXTRAN = 4)
C                                       Y WCS transformation
      INTEGER    GYTRAN
      PARAMETER (GYTRAN = 5)
C                                       Clip at viewport ?
      INTEGER    GRCLIP
      PARAMETER (GRCLIP = 6)
C                                       Draw title ?
      INTEGER    GLDTTL
      PARAMETER (GLDTTL = 34)
C                                       Viewport aspect ratio
      INTEGER    GASPCT
      PARAMETER (GASPCT = 38)
C                                       Character height in NDC units
      INTEGER    GCHRSZ
      PARAMETER (GCHRSZ = 39)
C                                       Polyline type
      INTEGER    GPLTYP
      PARAMETER (GPLTYP = 14)
C                                       Polymarker line type
      INTEGER    GPMTYP
      PARAMETER (GPMTYP = 18)
C                                       Text parameters:
C                                       ----------------
C                                       Text up direction
      INTEGER    GTXTUP
      PARAMETER (GTXTUP = 22)
C                                       Text size
      INTEGER    GTXSIZ
      PARAMETER (GTXSIZ = 23)
C                                       Text path
      INTEGER    GTXPTH
      PARAMETER (GTXPTH = 24)
C                                       Horizontal text justification
      INTEGER    GTHJST
      PARAMETER (GTHJST = 26)
C                                       Text vertical justification
      INTEGER    GTVJST
      PARAMETER (GTVJST = 27)
C                                       Text font
      INTEGER    GTFONT
      PARAMETER (GTFONT = 28)
C                                       Text quality
      INTEGER    GTQUAL
      PARAMETER (GTQUAL = 29)
C                                       X axis parameters:
C                                       ------------------
C                                       Draw x axes
      INTEGER    GLXDAX
      PARAMETER (GLXDAX = 101)
C                                       Set x axis position ?
      INTEGER    GLXSAP
      PARAMETER (GLXSAP = 102)
C                                       WCS of bottom x axis
      INTEGER    GLXAX1
      PARAMETER (GLXAX1 = 103)
C                                       WCS of top x axis
      INTEGER    GLXAX2
      PARAMETER (GLXAX2 = 104)
C                                       Draw x grid ?
      INTEGER    GLXDRG
      PARAMETER (GLXDRG = 105)
C                                       Round ends of x axis ?
      INTEGER    GLXRND
      PARAMETER (GLXRND = 106)
C                                       Label x axis ?
      INTEGER    GLXLBL
      PARAMETER (GLXLBL = 107)
C                                       Draw x axis ticks ?
      INTEGER    GLXDRT
      PARAMETER (GLXDRT = 109)
C                                       Label x axis ticks ?
      INTEGER    GLXLTK
      PARAMETER (GLXLTK = 110)
C                                       Number of major x axis ticks
      INTEGER    GLXNMJ
      PARAMETER (GLXNMJ = 111)
C                                       Number of minor x axis ticks
      INTEGER    GLXNMN
      PARAMETER (GLXNMN = 112)
C                                       Y axis parameters:
C                                       ------------------
C                                       Draw y axes
      INTEGER    GLYDAX
      PARAMETER (GLYDAX = 201)
C                                       Set y axis position ?
      INTEGER    GLYSAP
      PARAMETER (GLYSAP = 202)
C                                       WCS of bottom y axis
      INTEGER    GLYAX1
      PARAMETER (GLYAX1 = 203)
C                                       WCS of top y axis
      INTEGER    GLYAX2
      PARAMETER (GLYAX2 = 204)
C                                       Draw y grid ?
      INTEGER    GLYDRG
      PARAMETER (GLYDRG = 205)
C                                       Round ends of y axis ?
      INTEGER    GLYRND
      PARAMETER (GLYRND = 206)
C                                       Label y axis ?
      INTEGER    GLYLBL
      PARAMETER (GLYLBL = 207)
C                                       Draw y axis ticks ?
      INTEGER    GLYDRT
      PARAMETER (GLYDRT = 209)
C                                       Label y axis ticks ?
      INTEGER    GLYLTK
      PARAMETER (GLYLTK = 210)
C                                       Number of major y axis ticks
      INTEGER    GLYNMJ
      PARAMETER (GLYNMJ = 211)
C                                       Number of minor y axis ticks
      INTEGER    GLYNMN
      PARAMETER (GLYNMN = 212)
C                                       WCS transformations:
C                                       --------------------
C                                       Linear WCS transformation
      INTEGER    GWTLIN
      PARAMETER (GWTLIN = 0)
C                                       Log (10) WCS transformation
      INTEGER    GWTLOG
      PARAMETER (GWTLOG = 1)
C                                       Log (e) WCS transformation
      INTEGER    GWELOG
      PARAMETER (GWELOG = 2)
C                                       Polyline types:
C                                       ---------------
C                                       Clear (erase) polyline type
      INTEGER    GPLCLR
      PARAMETER (GPLCLR = 0)
C                                       Solid polyline type
      INTEGER    GPLSLD
      PARAMETER (GPLSLD = 1)
C                                       Dashed polyline type
      INTEGER    GPLDSH
      PARAMETER (GPLDSH = 2)
C                                       Dotted polyline type
      INTEGER    GPLDOT
      PARAMETER (GPLDOT = 3)
C                                       Dot-dashed polyline type
      INTEGER    GPLDTD
      PARAMETER (GPLDTD = 4)
C                                       Histogram styles:
C                                       -----------------
C                                       Step (pseudo-histogram) curve
      INTEGER    GHCSTP
      PARAMETER (GHCSTP = 2)
C                                       Bar graph (histogram)
      INTEGER    GHCBAR
      PARAMETER (GHCBAR = 3)
C                                       Polymarker types:
C                                       -----------------
C                                       Point marker type
      INTEGER    GMTPNT
      PARAMETER (GMTPNT = 0)
C                                       Box marker type
      INTEGER    GMTBOX
      PARAMETER (GMTBOX = 2)
C                                       Plus marker type
      INTEGER    GMTPLS
      PARAMETER (GMTPLS = 4)
C                                       Cross marker type
      INTEGER    GMTCRS
      PARAMETER (GMTCRS = 8)
C                                       Diamond marker type
      INTEGER    GMTDMD
      PARAMETER (GMTDMD = 16)
C                                       Horizontal line marker type
      INTEGER    GMTHLN
      PARAMETER (GMTHLN = 32)
C                                       Vertical line marker type
      INTEGER    GMTVLN
      PARAMETER (GMTVLN = 64)
C                                       Horizontal error bar marker type
      INTEGER    GMTHEB
      PARAMETER (GMTHEB = 128)
C                                       Vertical error bar marker type
      INTEGER    GMTVEB
      PARAMETER (GMTVEB = 256)
C                                       Circle marker type
      INTEGER    GMTCIR
      PARAMETER (GMTCIR = 512)
C                                       Text attributes:
C                                       ----------------
C                                       Center text
      INTEGER    GTACEN
      PARAMETER (GTACEN = 1)
C                                       Horizontally justify text left
      INTEGER    GTALFT
      PARAMETER (GTALFT = 2)
C                                       Horizontally justify text right
      INTEGER    GTARGT
      PARAMETER (GTARGT = 3)
C                                       Vertically justify text at top
      INTEGER    GTAVTP
      PARAMETER (GTAVTP = 6)
C                                       Vertically justify text at top
      INTEGER    GTAVBT
      PARAMETER (GTAVBT = 7)
C                                       Text path up
      INTEGER    GTAPUP
      PARAMETER (GTAPUP = 4)
C                                       Text path down
      INTEGER    GTAPDN
      PARAMETER (GTAPDN = 5)
C                                       Roman text font
      INTEGER    GTAROM
      PARAMETER (GTAROM = 8)
C                                       Greek text font
      INTEGER    GTAGRK
      PARAMETER (GTAGRK = 9)
C                                       Italic text font
      INTEGER    GTAITL
      PARAMETER (GTAITL = 10)
C                                       Bold text font
      INTEGER    GTABLD
      PARAMETER (GTABLD = 11)
C                                       Low quality (hardware) text
      INTEGER    GTALOW
      PARAMETER (GTALOW = 12)
C                                       High quality (software) text
      INTEGER    GTAHGH
      PARAMETER (GTAHGH = 14)
C                                       Graphics reset parameters:
C                                       --------------------------
C                                       Reset all parameters: 777 octal
      INTEGER    GRSETA
      PARAMETER (GRSETA = 511)
C                                       Reset GIO parameters
      INTEGER    GRSGIO
      PARAMETER (GRSGIO = 1)
C                                       Reset WCS parameters
      INTEGER    GRSWCS
      PARAMETER (GRSWCS = 2)
C                                       Reset GLABAX parameters
      INTEGER    GRSGLB
      PARAMETER (GRSGLB = 4)
C-----------------------------------------------------------------------------
C Status return codes that match in SPP and Fortran includes
C
C These codes all begin with the letters `US'.  The U identifies the parameter
C as being part of the F77/VOS interface, and the S signifies that these are
C status codes.  The third letter identifies the general category of function
C from which the code is returned.  These third characters are :
C
C     C  -  CL operations
C     G  -  graphics functions
C     H  -  header parameter I/O functions
C     I  -  image-related functions
C     M  -  dynamic memory
C     T  -  template processing
C     W  -  terminal/logfile I/O
C
C The first two parameters are special cases; the convention regarding the
C third character does not apply to these parameters.
C
C                                       End-of-file indicator
      INTEGER    USEOF
      PARAMETER (USEOF  =   -2)
C                                       Successful return
      INTEGER    USOK
      PARAMETER (USOK   =    0)
C                                       CL parameter access status messages:
C                                       ------------------------------------
C                                       End of file on CL string
      INTEGER    USCEOF
      PARAMETER (USCEOF =   -1)
C                                       CL end of list
      INTEGER    USCEOL
      PARAMETER (USCEOL =   -1)
C                                       CL parameter not found
      INTEGER    USCPNF
      PARAMETER (USCPNF =    1)
C                                       CL parameter has bad data type
      INTEGER    USCPDT
      PARAMETER (USCPDT =    2)
C                                       CL parameter is undefined
      INTEGER    USCPUD
      PARAMETER (USCPUD =    3)
C                                       Error putting parameter into CL
      INTEGER    USCPUT
      PARAMETER (USCPUT =    5)
C                                       First element is bad (must be > 0)
      INTEGER    USCBFE
      PARAMETER (USCBFE =    6)
C                                       Number of elements is bad (must be > 0)
      INTEGER    USCBNE
      PARAMETER (USCBNE =    7)
C                                       Image access status messages:
C                                       -----------------------------
C                                       Error opening existing file
      INTEGER    USIOPO
      PARAMETER (USIOPO =   10)
C                                       Error opening new file
      INTEGER    USIOPN
      PARAMETER (USIOPN =   11)
C                                       Error opening new copy file
      INTEGER    USIOPC
      PARAMETER (USIOPC =   12)
C                                       Invalid NAXIS parameter
      INTEGER    USINAX
      PARAMETER (USINAX =   13)
C                                       Invalid axis length parameter
      INTEGER    USIDIM
      PARAMETER (USIDIM =   14)
C                                       Invalid image pixel data type
      INTEGER    USIPDT
      PARAMETER (USIPDT =   15)
C                                       Error returning offset to image data
      INTEGER    USIOFS
      PARAMETER (USIOFS =   16)
C                                       Invalid image access mode
      INTEGER    USIACM
      PARAMETER (USIACM =   17)
C                                       Error closing existing file
      INTEGER    USICLO
      PARAMETER (USICLO =   18)
C                                       Not an image !
      INTEGER    USIBIM
      PARAMETER (USIBIM =   19)
C                                       Bad section specification for image
      INTEGER    USISEC
      PARAMETER (USISEC =   20)
C                                       Bad extension for image
      INTEGER    USIEXT
      PARAMETER (USIEXT =   21)
C                                       Illegal section on new image/new copy
      INTEGER    USIISC
      PARAMETER (USIISC =   22)
C                                       Error reading pixel file
      INTEGER    USIPRD
      PARAMETER (USIPRD =   23)
C                                       Error writing to pixel file
      INTEGER    USIPWR
      PARAMETER (USIPWR =   24)
C                                       Get unknown image parameter
      INTEGER    USIGUP
      PARAMETER (USIGUP =   25)
C                                       Set unknown image parameter
      INTEGER    USISUP
      PARAMETER (USISUP =   26)
C                                       Terminal I/O status messages:
C                                       -----------------------------
C                                       Error writing to STDOUT
      INTEGER    USWOWR
      PARAMETER (USWOWR =   30)
C                                       Error writing to the log file
      INTEGER    USWLWR
      PARAMETER (USWLWR =   31)
C                                       Invalid message destination
      INTEGER    USWBMD
      PARAMETER (USWBMD =   32)
C                                       Invalid message priority
      INTEGER    USWBMP
      PARAMETER (USWBMP =   33)
C                                       Header I/O status messages:
C                                       ---------------------------
C                                       Header parameter not found
      INTEGER    USHPNF
      PARAMETER (USHPNF =   40)
C                                       Illegal data type for header parameter
      INTEGER    USHPDT
      PARAMETER (USHPDT =   41)
C                                       Out of space in image header
      INTEGER    USHOOS
      PARAMETER (USHOOS =   42)
C                                       Illegal attempt to redefine parameter
      INTEGER    USHPDF
      PARAMETER (USHPDF =   43)
C                                       Bad input template (for parm tpt ??)
      INTEGER    USTBIT
      PARAMETER (USTBIT =   45)
C                                       Bad input first element
      INTEGER    USHBFE
      PARAMETER (USHBFE =   46)
C                                       Bad input number of elements
      INTEGER    USHBNE
      PARAMETER (USHBNE =   47)
C                                       Cannot delete standard keyword
      INTEGER    USHPDE
      PARAMETER (USHPDE =   48)
C                                       Cannot delete nonexistent keyword
      INTEGER    USHPNX
      PARAMETER (USHPNX =   49)
C                                       Image does not exist
      INTEGER    USINF
      PARAMETER (USINF  =   50)
C                                       Error deleting image
      INTEGER    USIDEL
      PARAMETER (USIDEL =   51)
C                                       Error renaming image
      INTEGER    USIREN
      PARAMETER (USIREN =   52)
C                                       Illegal image header
      INTEGER    USIBIH
      PARAMETER (USIBIH =   53)
C                                       Image descriptor mismatch
      INTEGER    USIDMM
      PARAMETER (USIDMM =   54)
C                                       Template processing status messages:
C                                       ------------------------------------
C                                       Error opening template
      INTEGER    USTOPE
      PARAMETER (USTOPE =   60)
C                                       Error getting next image from template
      INTEGER    USTGNX
      PARAMETER (USTGNX =   61)
C                                       Too many patterns in template filename
      INTEGER    USTMXP
      PARAMETER (USTMXP =   62)
C                                       Illegal filename template
      INTEGER    USTBFT
      PARAMETER (USTBFT =   63)
C                                       Graphics processing status messages:
C                                       ------------------------------------
C                                       Bad graphics device
      INTEGER    USGDEV
      PARAMETER (USGDEV =   70)
C                                       Bad graphics mode
      INTEGER    USGMOD
      PARAMETER (USGMOD =   71)
C                                       Error opening graphics device
      INTEGER    USGOPE
      PARAMETER (USGOPE =   72)
C                                       Error closing graphics device
      INTEGER    USGCLO
      PARAMETER (USGCLO =   73)
C                                       Error clearing graphics device
      INTEGER    USGCLR
      PARAMETER (USGCLR =   74)
C                                       Error setting window
      INTEGER    USGSWN
      PARAMETER (USGSWN =   75)
C                                       Error getting window
      INTEGER    USGGWN
      PARAMETER (USGGWN =   76)
C                                       Error plotting a polyline
      INTEGER    USGPLN
      PARAMETER (USGPLN =   77)
C                                       Error plotting a sequence of markers
      INTEGER    USGMRK
      PARAMETER (USGMRK =   78)
C                                       Invalid number of points
      INTEGER    USGBNP
      PARAMETER (USGBNP =   79)
C                                       Bad marker (size is incorrect)
      INTEGER    USGBMK
      PARAMETER (USGBMK =   80)
C                                       Bad axes specification
      INTEGER    USGBAX
      PARAMETER (USGBAX =   81)
C                                       Bad justification for test
      INTEGER    USGBJS
      PARAMETER (USGBJS =   82)
C                                       Bad character size
      INTEGER    USGBCS
      PARAMETER (USGBCS =   83)
C                                       Bad orientation (for text)
      INTEGER    USGBOR
      PARAMETER (USGBOR =   84)
C                                       Error writing label
      INTEGER    USGLAB
      PARAMETER (USGLAB =   85)
C                                       Error writing text
      INTEGER    USGTXT
      PARAMETER (USGTXT =   86)
C                                       Error setting viewport
      INTEGER    USGSVP
      PARAMETER (USGSVP =   87)
C                                       Error getting viewport
      INTEGER    USGGVP
      PARAMETER (USGGVP =   88)
C                                       Error setting GIO parameter
      INTEGER    USGPP
      PARAMETER (USGPP  =   89)
C                                       Error getting GIO parameter
      INTEGER    USGGP
      PARAMETER (USGGP  =   90)
C                                       Bad reset code (in reset function)
      INTEGER    USGBRC
      PARAMETER (USGBRC =   91)
C                                       Error writing cursor
      INTEGER    USGCWR
      PARAMETER (USGCWR =   92)
C                                       Error in high level plot routine
      INTEGER    USGPLT
      PARAMETER (USGPLT =   93)
C                                       Error in high level error bar routine
      INTEGER    USGERB
      PARAMETER (USGERB =   94)
C                                       Error in high level histogram plot
      INTEGER    USGHGM
      PARAMETER (USGHGM =   95)
C                                       Bad cursor call
      INTEGER    USGBCC
      PARAMETER (USGBCC =   96)
C                                       Error in scale procedure
      INTEGER    USGSCA
      PARAMETER (USGSCA =   97)
C                                       Dynamic memory status messages:
C                                       -------------------------------
C                                       Invalid dynamic memory data type
      INTEGER    USMDT
      PARAMETER (USMDT  =  100)
C                                       Error allocating dynamic memory
      INTEGER    USMALL
      PARAMETER (USMALL =  101)
C                                       Error freeing dynamic memory
      INTEGER    USMFRE
      PARAMETER (USMFRE =  102)
C
C   Utility routines symbolic codes
C
C   The following symbolic codes are use in
C   routine 'UUDTOC'. (Convert a double precision
C   number into a string with a given format).
C   The symbolic format are as follows:


C                                       Decimal representaion ('d')
      INTEGER    FMTDEC
      PARAMETER (FMTDEC =  100)
C                                       Exponential form ('e')
      INTEGER    FMTEXP
      PARAMETER (FMTEXP =  101)
C                                       Fixed floating point ('f')
      INTEGER    FMTFIX
      PARAMETER (FMTFIX =  102)
C                                       General floating point(E,F) ('g')
      INTEGER    FMTGEN
      PARAMETER (FMTGEN =  103)
C                                       Hour:minutes:seconds representation('h')
      INTEGER    FMTHMS
      PARAMETER (FMTHMS =  104)
C                                       Minutes:seconds representation ('m')
      INTEGER    FMTMIS
      PARAMETER (FMTMIS =  109)
C                                       Octal representation ('o')
      INTEGER    FMTOCT
      PARAMETER (FMTOCT =  111)
C                                       Hex representation ('x')
      INTEGER    FMTHEX
      PARAMETER (FMTHEX =  120)
C
C  Error codes from utility routines
C
C     					Not a number
      INTEGER    USUNAN
      PARAMETER  (USUNAN =  110)
C     					Overflow
      INTEGER    USUOVF
      PARAMETER  (USUOVF =  111)
C     					String not long enough to contain
C					output character decimal.
      INTEGER    USUNSP
      PARAMETER  (USUNSP =  112)
C-------------------------------------- End of IRAF77.INC
