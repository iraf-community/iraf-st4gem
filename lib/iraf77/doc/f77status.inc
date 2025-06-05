C IRAF77 -- Symbolic parameter definitions for F77/VOS interface
C
C
C Status return codes that match in SPP and Fortran includes
C
C                                       End-of-file indicator; same as IRAF
      INTEGER    USEOF
      PARAMETER (USEOF  =   -2)
C                                       Successful return
      INTEGER    USOK
      PARAMETER (USOK   =    0)
C
C CL parameter access status messages : codes 1-9
C
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
C
C Image access status messages : codes 10-29
C
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
C
C Terminal I/O status messages : codes 30-39
C
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
C
C Header I/O status messages : 40-59
C
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
C
C Template processing status messages : codes 60-69
C
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
C                                       Error closing template
      INTEGER    USTCLS
      PARAMETER (USTCLS =   64)
C
C Graphics processing status messages : codes 70-99
C
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
C
C Dynamic memory status messages : codes 100-109
C
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
C Utility character/number conversion status messages : 110-119
C
C                                       Not a number
      INTEGER    USUNAN
      PARAMETER (USUNAN =  110)
C                                       Overflow
      INTEGER    USUOVF
      PARAMETER (USUOVF =  111)
C                                       Field width not wide enough
      INTEGER    USUNSP
      PARAMETER (USUNSP =  112)
