        SUBROUTINE RDWAVE (

*  inputs:

     :                     VERB,
     :                     WAVTBL,

*  outputs:

     :                     NWAVE,
     :                     WAVE
     :                    )
*
*  Keyphrase:
*  ----------
*  Read the requested wavelength set from the wavelength table
*
*  Description:
*  ------------
*  This subroutine opens the wavelength table, WAVTBL, and reads the
*  requested wavelength set into a dynamically allocated array. Null
*  values are removed from the array, which is returned to the calling
*  routine.
*
*  FORTRAN Name: RDWAVE.FOR
*
*  Keywords of Accessed Files and Tables:
*  --------------------------------------
*  Name              I/O       Description/comment
*  'WAVTABLE'         I        SDAS table contains requested wavelength set
*
*  Subroutines Called:
*  -------------------
*  CDBS:
*       ERRMSG, DTLMSG 
*  SDAS:
*       UTTOPN, UTPGTI, UTCFND, UTCGTR, UDMGET, UTTCLO
*  Others: None
*
*  History:
*  --------
*  Version     Date          Author        Description
*     1      04-02-87        K.-T. Yung    code
*     2      06-02-88        B.Simon       revised for modularity
*
*-----------------------------------------------------------------------------
*
*==: IRAF.H equivalents
*
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
      INTEGER          MEMI(1)
      REAL             MEMR(1)
      DOUBLE PRECISION MEMD(1)
      COMPLEX          MEMX(1)
      EQUIVALENCE (MEMB, MEMI, MEMR, MEMD, MEMX)
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
C-------------------------------------- End of IRAF77.INC
*
*== String array dimensions:
*
*                               --maximum number of components
        INTEGER MAXCMP,
*                               --maximum number of keywords
     :          MAXKEY,
*                               --maximum component per mode
     :          MAXLST
*
        PARAMETER ( MAXCMP = 2000 )
*
        PARAMETER ( MAXKEY = 100 )
*
        PARAMETER ( MAXLST = 100 )
*
*== Miscellaneous parameters:
*
*                               -- success indicator
        INTEGER         OK,
*                               -- error indicator
     :                  NOTOK
*
        PARAMETER ( OK = 0 )
*
        PARAMETER ( NOTOK = 1 )
*
*== input: 
*                               --detail message switch
        LOGICAL         VERB
*                               --table contains requested wavelength set
        CHARACTER*64    WAVTBL
*
*== output: 
*                               --number of wavelengths in table
        INTEGER         NWAVE 
*                               --pointer to wavelength set
        INTEGER         WAVE
*
*== local:              
*                               --wavelength table pointer
        INTEGER         TPWAV, 
*                               --column indicator
     :                  CIWAV, 
*                               --number of non-null wavelengths
     :                  MWAVE,
*                               --SDAS return status
     :                  STATUS
*                               --array of flags (true if undefined)
        INTEGER         NULFLG
*
*---------------------------------------------------------------------------

*       Get number of rows in wavelength set table

        CALL UTTOPN( WAVTBL, RDONLY, TPWAV, STATUS )
        CALL ERRMSG( STATUS, WAVTBL, 'Unable to open table @' )
        CALL DTLMSG( VERB, WAVTBL, 'Wavelength table @ opened' )

        CALL UTPGTI( TPWAV, TBNROW, NWAVE, STATUS )
        CALL ERRMSG( STATUS, WAVTBL,
     :               'Unable to get number of rows in @' )

*       Allocate the dynamic memory to hold the wavelength table column

        CALL UDMGET ( NWAVE, TYBOOL, NULFLG, STATUS )
        CALL UDMGET ( NWAVE, TYREAL, WAVE, STATUS )
        CALL ERRMSG ( STATUS, ' ', 'Unable to allocate dynamic memory' )

*       Get the requested wavelengths from WAVTBL, place them in array WAVE(I)

        CALL UTCFND( TPWAV, 'WAVELENGTH', 1, CIWAV , STATUS )
        CALL ERRMSG( STATUS, WAVTBL,
     :               'Column WAVELENGTH not found in @' )

        CALL UTCGTR( TPWAV, CIWAV, 1, NWAVE, MEMR(WAVE), MEMB(NULFLG),
     :               STATUS )
        CALL ERRMSG( STATUS, WAVTBL,
     :               'Unable to get requested wavelengths from @' )

*       Move all null elements to the end of the array

        CALL DNULLR( MEMB(NULFLG), NWAVE, MWAVE, MEMR(WAVE) )
        IF ( MWAVE .EQ. 0) THEN
            CALL ERRMSG( NOTOK, WAVTBL,
     :                  'Wavelength values are all undefined in @' )
        ELSE
            NWAVE = MWAVE
        END IF

        CALL UTTCLO( TPWAV, STATUS )
        CALL ERRMSG( STATUS, WAVTBL, 'Unable to close table @' )

        RETURN
        END
 
                                                               
                                                               
                                                               
                                                               
                                                               
                                                               
                                                              

 
