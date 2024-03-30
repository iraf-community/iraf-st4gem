        SUBROUTINE PUTBND (VERB, NOWAVE, MODSTR, NWAVE, WAVE, 
     :		           GRFTBL, CMPTBL, FLTTAB, STATUS)
*  Worthless comment to see what happens
*
*  Keyphrase:
*  ----------
*  Calculate the instrument passband for all components in a graph path 
*
*  Description:
*  ------------
*  This subroutine first determines the list of optical component mnemonics
*  from the mode string and the graph table by calling PATHF.  Then it gets
*  the corresponding file names from the component table by calling CMPFIL.
*  The wavelengths at which the throughput is desired are pointed to by WAVE,
*  i.e. MEMR(WAVE) through MEMR(WAVE + NWAVE-1).
*  This subroutine then computes the instrument passband and error on this
*  wavelength set by calling THUSET.  Finally, it writes the wavelength set,
*  throughput, and error to the instrument passband table, FLTTAB.
*
*  FORTRAN Name: PUTBND.FOR
*
*  Keywords of Accessed Files and Tables:
*  --------------------------------------
*  NAME              I/O       Description/comment
*  GRFTBL             I        HST instrument graph table
*  CMPTBL             I        HST component name table
*  FLTTAB             O        Output table with instrument passband
*
*  Subroutines Called:
*  -------------------
*  CDBS:
*       ERRMSG, PATHF, CMPFIL, WAVSET, THUSET, FUNIT, PIVLAM,
*       PHOZPT, RMSLAM
*  SDAS:
*       UCLGST, UCLGSB, UTTINN, UTCDEF, UTTCRE, UTCPTR,
*       UTHADR, UTHADT, UTTCLO, UDMGET, UDMFRE
*  OTHERS: NONE
*
*  History:
*  --------
*  Version     Date          Author        Description
*     1      08-16-89        B.Simon       new synphot library interface
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
C                                       File access modes, etc.:
C                                       ------------------------
C                                       Read only
      INTEGER    RDONLY
      PARAMETER (RDONLY = 1)
C                                       Read/Write
      INTEGER    RDWRIT
      PARAMETER (RDWRIT = 2)
C                                       New file
      INTEGER    NEWFIL
      PARAMETER (NEWFIL = 5)
C                                       New copy image
      INTEGER    NEWCPY
      PARAMETER (NEWCPY = 7)
C                                       Single-precision real
      INTEGER    TYREAL
      PARAMETER (TYREAL = 6)
C                                       Double-precision real
      INTEGER    TYDOUB
      PARAMETER (TYDOUB = 7)
C------------------------------------------------------------------------------
C Definitions for parameters relevant to table I/O.
C
C                                       Maximum number of header parameters
      INTEGER    TBMXPR
      PARAMETER (TBMXPR = 6)
C                                       Number of rows written to
      INTEGER    TBNROW
      PARAMETER (TBNROW = 21)
C                                       Number of columns defined
      INTEGER    TBNCOL
      PARAMETER (TBNCOL = 22)
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
*== input:
*                               --verbose message flag
        LOGICAL         VERB,
*                               --no wavelength array flag
     :                  NOWAVE
*                               --string with mode keywords or mode filename
        CHARACTER*(*)   MODSTR,
*                               --name of output passband table
     :                  FLTTAB
*                               --number of wavelengths
        INTEGER         NWAVE

*
*== input/output:
*                               --pointer to (real) array of wavelengths
        INTEGER         WAVE
*                               --name of instrument graph table
        CHARACTER*(*)   GRFTBL,
*                               --name of component list lookup table
     :                  CMPTBL

*
*== output:
*                               --sdas return status
        INTEGER         STATUS
*
*== local:              
*                               --array contains file names of thruput 
        CHARACTER*64    THUTBL(MAXLST)
*                               --component mnemonics of graph links
        CHARACTER*20    LIST(MAXLST)
*                               --instrument passband table pointer

        INTEGER         TPFLT,
*                               --column indicators
     :                  CIFLT(3),
*                               --total number of components present
     :                  NLIST,
*                               --number of components with throughput tables
     :                  NTHU,
*                               --old value of status
     :                  OSTAT
*                               --passband values at wavelength

        INTEGER         FILT,
*                               --passband errors at wavelength
     :                  FLTERR

*                               --photometric parameters
        REAL            PHOT(5)

*                               --column names of instrument passband table
        CHARACTER*20    FLTCOL(3),
*                               --units of output instrument passband table
     :                  UNITS(3),
*                               --formats of passband table columns
     :                  FORMAT(3)

*                               --data type of passband table columns
        INTEGER         DTTYPE(3)
*
        DATA    FLTCOL  /'WAVELENGTH', 'THROUGHPUT', 'ERROR' /
*
        DATA    UNITS   /'ANGSTROMS', ' ', ' ' /
*
        DATA    FORMAT  /'F12.3', 'F12.7', 'F12.7' /
*
        DATA    DTTYPE  /TYREAL, TYREAL, TYREAL /
*
*==: functions:
*       
*                               --calculate the unit response
        REAL            FUNIT,
*                               --calculate the pivot wavelength
     :                  PIVLAM,
*                               --calculate the equivalent RMS wavelength
     :                  RMSLAM,
*                               --get a synthetic photometry constant
     :                  SYNVAL
*
*---------------------------------------------------------------------------

*       Get component keywords corresponding to observation mode

        CALL PATHF (VERB, MODSTR, MAXLST, GRFTBL, NLIST, LIST)

*       Get the desired optical mnemonic names from CMPTBL

        CALL CMPFIL (VERB, NLIST, LIST, CMPTBL, NTHU, THUTBL)

*       Create the wavelength set at which the passband functions will
*       be evaluated if it was not supplied by the user

        IF (NOWAVE) THEN
            CALL WAVSET (NTHU, THUTBL, NWAVE, WAVE)
	ENDIF

*       Allocate dynamic memory for throughput arrays

        CALL UDMGET (NWAVE, TYREAL, FILT, STATUS)
        CALL UDMGET (NWAVE, TYREAL, FLTERR, STATUS)
        CALL ERRMSG (STATUS, ' ', 'Unable to allocate virtual memory')

*       Compute the combined throughput

        CALL THUSET (NTHU, THUTBL, NWAVE, WAVE, FILT, FLTERR)

*       Evaluate passband functions of components on wavelength set

        PHOT(1) = FUNIT (NWAVE, MEMR(WAVE), MEMR(FILT), 'flam')
        PHOT(2) = SYNVAL ('stmag')
        PHOT(3) = PIVLAM (NWAVE, MEMR(WAVE), MEMR(FILT))
        PHOT(4) = RMSLAM (NWAVE, MEMR(WAVE), MEMR(FILT))
        PHOT(5) = SYNVAL ('hstarea')

        CALL UTTINN (FLTTAB, TPFLT, STATUS)
        IF (STATUS .NE. OK) GOTO 900

        CALL UTPPTI (TPFLT, TBMXPR, 10, STATUS)
        IF (STATUS .NE. OK) GOTO 900

        CALL UTCDEF (TPFLT, FLTCOL, UNITS, FORMAT, DTTYPE, 3,
     :               CIFLT, STATUS)
        IF (STATUS .NE. OK) GOTO 900

        CALL UTTCRE (TPFLT, STATUS)
        IF (STATUS .NE. OK) GOTO 900

*       Put computed values in table header.

        CALL UTHADR (TPFLT, 'photflam', PHOT(1), STATUS)
        CALL UTHADR (TPFLT, 'photzpt', PHOT(2), STATUS)
        CALL UTHADR (TPFLT, 'photplam', PHOT(3), STATUS)
        CALL UTHADR (TPFLT, 'photbw', PHOT(4), STATUS)
        CALL UTHADR (TPFLT, 'aperarea', PHOT(5), STATUS)
        IF (STATUS .NE. OK) GOTO 900

        CALL UTHADT (TPFLT, 'obsmode', MODSTR, STATUS)
        CALL UTHADT (TPFLT, 'grftable', GRFTBL, STATUS)
        CALL UTHADT (TPFLT, 'cmptable', CMPTBL, STATUS)
        IF (STATUS .NE. OK) GOTO 900

        CALL UTCPTR (TPFLT, CIFLT(1), 1, NWAVE, MEMR(WAVE), STATUS)
        CALL UTCPTR (TPFLT, CIFLT(2), 1, NWAVE, MEMR(FILT), STATUS)
        CALL UTCPTR (TPFLT, CIFLT(3), 1, NWAVE, MEMR(FLTERR), STATUS)
        IF (STATUS .NE. OK) GOTO 900

        CALL DTLMSG (VERB, FLTTAB, 'Instrument passband is in table @')
        CALL UTTCLO (TPFLT, STATUS)

 900    IF (STATUS .EQ. OK) THEN
	    OSTAT = OK
        ELSE
	    OSTAT = NOTOK
        ENDIF

*       Free dynamic memory 

        CALL UDMFRE (FILT, TYREAL, STATUS)
        CALL UDMFRE (FLTERR, TYREAL, STATUS)
        CALL ERRMSG (STATUS, ' ', 'Unable to free virtual memory')

        STATUS = OSTAT

        RETURN
        END
