C     Copyright restrictions apply - see st4gem$copyright.stsdas 
C 
      SUBROUTINE PLUCY   !I
CM    PROGRAM MPLUCY      !M
C++
C
C PLUCY.F Version 0.7
C
C Lucy point-sources+smooth extended objects restoration.
C Experimental version for many PSFs, acceleration and a floating prior.
C
C Implementation Details:
C
C This version uses the F77/VOS interface to IRAF parameters and
C data structures.
C
C Apart from the F77/VOS standard libraries the only other external
C dependencies are:
C
C fftetc.f - a set of FFT support routines including several
C            FFT codes.
C
C x_plucy.x - the one line of SPP needed for linking purposes.
C
C timtem.x - routines needed for image list handling from F77/VOS
C
C Author:
C
C Richard Hook, ST-ECF, e-mail: rhook@eso.org (Internet)
C (Based on ideas and code by Leon Lucy)
C
C History:
C
C First attempts, Jan 1993
C Added renormalisation and non-negativity checking - July 1993
C Started multi PSF version - Aug 93
C Floating prior added - Aug 93
C Acceleration added - Aug 93
C Acceleration for stationary prior and other additions - Sep/Oct 93
C Tests of automatically corrected PSF - Oct 93
C Bugfix for version 0.5D - Mar 94
C Working on PSF determination enhancements - Aug 94
C Version 0.6D, including regularization of PSF - Aug 94
C  This version also modifies and simplifies the normalisations
C Version 0.65 - bugfix for DIMS array - July 95
C
C Version 0.7 - removed options for varying PSF, removed problems
C               with character strings extending over more than one
C               line and a few other minor changes.
C
C Note:
C
C This version may also be compiled and linked to run under MIDAS
C as well as IRAF. Essentially all that is required is for the
C lines ending with !I to be commented out and the comment
C characters removed from lines starting CM. The result should
C then be linked with the additional library 'midint.f' which
C emulates (partially) the F77/VOS routines using MIDAS.
C
C--
      IMPLICIT NONE
                   
C Iraf global storage
      DOUBLE PRECISION MEMD(1) !I
      COMMON /MEM/MEMD         !I
                                 
C MIDAS global storage
CM    INTEGER MEMD(1)          !M
CM    COMMON /VMR/MEMD         !M
                                 
C Parameters
      INTEGER MAXSTR
      PARAMETER (MAXSTR=100000)
      INTEGER MAXPSF
      PARAMETER (MAXPSF=100)
      INTEGER USEOF
      PARAMETER (USEOF=-2)

C Local variables
      INTEGER DATTYP,NDIMS,NPSF
      INTEGER DATPNT,PSFPNT(MAXPSF)
      INTEGER ISTAT,NX,NY,N,NITER,DIMS(7)
      INTEGER PSFFFT(MAXPSF),WORK,REST,PS,PSM(MAXPSF)
      INTEGER W,CH,NNEG,CF,PHS,RP,CFP,PSFP,FPFFT,PHP
      INTEGER DPHP,DPHS
      INTEGER IDD,IDP(MAXPSF),IDDC,IDPS
      DOUBLE PRECISION T,PT,SCRATCH(40960),CS,SMB,TST
      DOUBLE PRECISION RMSRES,RESMAX,TT,TEM,MAGZERO
      DOUBLE PRECISION GW,HH,SS,STARS(5,MAXSTR),QQ
      INTEGER IRMAX,JRMAX,NSTARS
      INTEGER I,PSFLD,IMSTAT
      DOUBLE PRECISION ST,XL,VMIN,LV
      CHARACTER*80 DATA,PSFS,PSF(MAXPSF),FILE,OUTTAB
      CHARACTER*80 DECON,BACK
      CHARACTER*80 CHARS
      INTEGER NBX,NBY
      LOGICAL VERBOSE
      LOGICAL SUBSAM
      LOGICAL FPRIOR
      LOGICAL ACCEL 
C++
   
C Start of code
               
CM    CALL STSPRO('MPLUCY') !M
                              
C First get the boolean control logicals
      CALL UCLGSB('verbose',VERBOSE,ISTAT)
      CALL UCLGSB('accel',ACCEL,ISTAT)

      FPRIOR=.TRUE.

C Initialise the acceleration factor to one
      XL=1.0D0

C Initialise the maximised quantity (obj func.) to something v. small
      LV=-1.0E20

C Announce the version
      IF(VERBOSE) CALL UMSPUT('+ PLUCY Version 0.7 (March 98)', !I
     :            1,0,ISTAT)                                      !I
CM    IF(VERBOSE) CALL UMSPUT('+ MPLUCY Version 0.7 (March 98)', !M
CM   :            1,0,ISTAT)                                      !M
                            
C Get the magnitude zero point
      CALL UCLGSD('magzero',MAGZERO,ISTAT)
                                    
C Get the name of the input data image
      CALL UCLGST('data',DATA,ISTAT)
      CALL UIMOPN(DATA,1,IDD,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT('! Unable to open data image',
     :                    1,0,ISTAT)
         GO TO 99
      ELSE
         IF(VERBOSE) CALL UMSPUT('-Opening data file: '//DATA,
     :                   1,0,ISTAT)
      ENDIF
         
C Get the names of the PSF images and count them
C if the PSF is varying this is the input, initial PSF
C which need not be the same as the default
      CALL UCLGST('psfs',PSFS,ISTAT)
      CALL TIMOTP(PSFS,PSFLD,ISTAT)

      NPSF=1
      IMSTAT=0
      DO I=1,MAXPSF+1
         CALL TIMXTP(PSFLD,PSF(NPSF),IMSTAT)
         IF(IMSTAT.EQ.USEOF) GO TO 78
         CALL UIMOPN(PSF(NPSF),1,IDP(NPSF),ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT('! Unable to open PSF image',
     :                    1,0,ISTAT)
            GO TO 99
         ELSE
             IF(VERBOSE) 
     :          CALL UMSPUT('-Opening PSF file: '//PSF(NPSF),
     :                   1,0,ISTAT)
         ENDIF
         NPSF=NPSF+1
      ENDDO

 78   CONTINUE
      NPSF=NPSF-1

      CALL TIMCTP(PSFLD,ISTAT)

C Check for floating prior and multiple PSFs - this is not currently
C supported
      IF(NPSF.GT.1 .AND. (FPRIOR.OR.ACCEL)) THEN
         CALL UMSPUT('! Sorry, multiple PSFs cannot be used'//
     :               ' with a floating prior or acceleration.',
     :               1,0,ISTAT)
         GO TO 99
      ENDIF

C Get the fractional contribution of the entropy term
      CALL UCLGSD('fracent',CS,ISTAT)

C Do some checks on PSF size etc
C Get the shapes and sizes of the images
      CALL UIMGID(IDD,DATTYP,NDIMS,DIMS,ISTAT)
                                                    
      IF(NDIMS.EQ.1) THEN
         DIMS(2)=1
      ENDIF
           
      NX=DIMS(1)
      NY=DIMS(2)
                
C Check that the dimensions are EVEN, if not abort
      IF(MOD(NX,2).NE.0 .OR.
     :   (MOD(NY,2).NE.0 .AND. NY.NE.1)) THEN
         CALL UMSPUT('! Image dimensions must be even',
     :               1,0,ISTAT)
         GO TO 99
      ENDIF

C Same for the PSFs
      DO N=1,NPSF
         CALL UIMGID(IDP(N),DATTYP,NDIMS,DIMS,ISTAT)
                                                      
         IF(NDIMS.EQ.1) THEN
            DIMS(2)=1
         ENDIF
              
C Check they are OK
         IF(NX.NE.DIMS(1) .OR.
     :        NY.NE.DIMS(2)) THEN
            CALL UMSPUT('! Data/PSF size/shape error',1,0,ISTAT)
            WRITE(CHARS,'(''Data array size: '',2I5,'//
     :           ''', PSF['',I4,''] array size: '',2I5)')
     :           NX,NY,N,DIMS(1),DIMS(2)
            CALL UMSPUT(CHARS,1,0,ISTAT)
            GO TO 99
         ENDIF
      ENDDO          

C If we have a floating prior get the smoothing kernel sigma
C and create an image for it. Then take the FFT and remove the
C original image
      IF(FPRIOR) THEN
         CALL UCLGSD('skernel',GW,ISTAT)

         CALL UDMGET(NX*NY,7,PSFP,ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :           '! Unable to allocate memory for PSFP array',
     :                    1,0,ISTAT)
            GO TO 99
         ENDIF

         CALL UDMGET(NX*NY*2,7,FPFFT,ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :           '! Unable to allocate memory for FPFFT array',
     :                    1,0,ISTAT)
            GO TO 99
         ENDIF

C Put a central, circular, normalised Gaussian into the image
         CALL FGAUSS(MEMD(PSFP),NX,NY,GW)

         CALL DFILL(MEMD(PSFP),NX,NY,MEMD(FPFFT))

         CALL DFOURT(MEMD(FPFFT),DIMS,2,-1,0,SCRATCH)

         CALL UDMFRE(PSFP,7,ISTAT)

C We also get space for the other images we will need
         CALL UDMGET(NX*NY,7,RP,ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :           '! Unable to allocate memory for RP array',
     :                    1,0,ISTAT)
           GO TO 99
         ENDIF

         CALL UDMGET(NX*NY,7,CFP,ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :           '! Unable to allocate memory for CFP array',
     :                    1,0,ISTAT)
           GO TO 99
         ENDIF

         CALL UDMGET(NX*NY,7,PHP,ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :           '! Unable to allocate memory for PHP array',
     :                    1,0,ISTAT)
           GO TO 99
         ENDIF
      ENDIF

C Even more arrays for the accelerated case
      IF(ACCEL) THEN

C This one is only needed for the floating prior case
         IF(FPRIOR) THEN
            CALL UDMGET(NX*NY,7,DPHP,ISTAT)
            IF(ISTAT.NE.0) THEN
               CALL UMSPUT(
     :         '! Unable to allocate memory for DPHP array',
     :                       1,0,ISTAT)
               GO TO 99
            ENDIF
         ENDIF

         CALL UDMGET(NX*NY,7,DPHS,ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :           '! Unable to allocate memory for DPHS array',
     :                    1,0,ISTAT)
            GO TO 99
         ENDIF
      ENDIF

C Get a scratch array
      CALL UDMGET(NX*NY,7,W,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :      '! Unable to allocate memory for working array',
     :            1,0,ISTAT)
         GO TO 99
      ENDIF
      CALL FILCON(MEMD(W),NX,NY,0.0D0)

C Get the name of the stars list
      CALL UCLGST('starlist',FILE,ISTAT)

      IF(FILE.NE.' ') THEN
         CALL GETSTR(FILE,STARS,NX,NY,NPSF,MEMD(W),
     :               MAXSTR,MAGZERO,
     :               NSTARS,ISTAT)
         IF(ISTAT.NE.0) GO TO 99

C Report the file name of the number of stars
         IF(VERBOSE) THEN
            WRITE(CHARS,'(''-Read star list file: '','//
     :                    'A20,'' ('',I5,'' stars).'')')
     :      FILE,NSTARS
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

C Also get the name for the output text file
         CALL UCLGST('outtab',OUTTAB,ISTAT)

         IF(OUTTAB.NE.' ') THEN
            OPEN(8,FILE=OUTTAB,STATUS='NEW',IOSTAT=ISTAT)
            IF(ISTAT.NE.0) THEN
               CALL UMSPUT(
     :      '! Unable to open output text list file',
     :                     1,0,ISTAT)
               GO TO 99
            ENDIF
         ENDIF
      ELSE
         NSTARS=0
         NPSF=1
      ENDIF

C Allocate space for the data
      CALL UDMGET(NX*NY,7,DATPNT,ISTAT)
      IF(ISTAT.NE.0) THEN
        CALL UMSPUT(
     :        '! Unable to allocate memory for data array',
     :                 1,0,ISTAT)
        GO TO 99
      ENDIF
             
C Read the data into the memory
      CALL UIGS2D(IDD,1,NX,1,NY,MEMD(DATPNT),ISTAT)
      IF(ISTAT.NE.0) THEN
        CALL UMSPUT(
     :        '! Unable to read in data input image',
     :                 1,0,ISTAT)
        GO TO 99
      ENDIF
             
C Check that the data is non-negative
      CALL ZAPNEG(MEMD(DATPNT),NX,NY,NNEG,VMIN,0.0D0)
      IF(NNEG.GT.0) THEN
          WRITE(CHARS,
     :    '(''! Warning, data image'','//
     :    ''' contains'',I7,'' negative values, these'','//
     :    ''' have been zeroed.'')') NNEG
         CALL UMSPUT(CHARS,1,0,ISTAT)
         WRITE(CHARS,'(''--Minimum value: '',d20.10)') VMIN
         CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF
             
C Allocate space for the PSFs
      DO N=1,NPSF
         CALL UDMGET(NX*NY,7,PSFPNT(N),ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :     '! Unable to allocate memory for PSF array',
     :              1,0,ISTAT)
           GO TO 99
         ENDIF
             
C Read the PSF into the memory
         CALL UIGS2D(IDP(N),1,NX,1,NY,MEMD(PSFPNT(N)),ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :      '! Unable to read in PSF input image',
     :               1,0,ISTAT)
            GO TO 99
         ENDIF
             
C Get working space arrays for the FFT
         CALL UDMGET(NX*NY*2,7,PSFFFT(N),ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :      '! Unable to allocate memory for PSF FFT array',
     :               1,0,ISTAT)
            GO TO 99
         ENDIF
      ENDDO
           
C Create the output data arrays
      CALL UCLGST('decon',DECON,ISTAT)
      CALL UIMCRE(DECON,7,2,DIMS,IDDC,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :     '! Unable to open output deconvolved image',
     :                 1,0,ISTAT)
         GO TO 99
      ENDIF
           
      CALL UCLGST('back',BACK,ISTAT)

      IF(BACK.NE.' ') THEN
         CALL UIMCRE(BACK,7,2,DIMS,IDPS,ISTAT)

         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :     '! Unable to open output background image',
     :                 1,0,ISTAT)
            GO TO 99
         ENDIF
      ENDIF
           
89    CONTINUE
              
C Get the sub-sampling values in X and Y and check that
C the X and Y dimensions are multiples of them.
      CALL UCLGSI('xsubsam',NBX,ISTAT)
      CALL UCLGSI('ysubsam',NBY,ISTAT)
                                      
      IF(MOD(NX,NBX).NE.0 .OR.
     :   MOD(NY,NBY).NE.0) THEN
         CALL UMSPUT('! Warning, arrays are not multiples'//
     :               ' of the sub-sampling',1,0,ISTAT)
      ENDIF
           
      IF(NBX.NE.1 .OR. NBY.NE.1) THEN
         SUBSAM=.TRUE.
      ELSE
         SUBSAM=.FALSE.
      ENDIF
           
C Check the sub-stepping is actually what is found in the data
      IF(SUBSAM) THEN
         CALL CHKSUB(MEMD(DATPNT),NX,NY,NBX,NBY,ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :    '! Sampling of data wrong - check x/ysubsam',
     :      1,0,ISTAT)
            GO TO 99
         ENDIF
      ENDIF

C Add up the total in the data image
      CALL TOTAL(MEMD(DATPNT),NX,NY,T)

C Normalise the image - this is necessary as the algorithm 
C assumes it has been done
      CALL MULC(MEMD(DATPNT),NX,NY,1.0D0/T,MEMD(DATPNT))

C Prepare the arrays for the FFT routine
C and do the FFTs of the PSF and rotated PSF
C
C At this point we also check that the PSFs are normalised to
C a total of 1 - if they aren't we warn the user but continue
                
C First check for negative values, if there are any set them
C to zero
      DO N=1,NPSF
         CALL ZAPNEG(MEMD(PSFPNT(N)),NX,NY,NNEG,VMIN,0.0D0)
         IF(NNEG.GT.0) THEN
            WRITE(CHARS,
     :       '(''! Warning, PSF '','//
     :       ''' contained '',I7,'' negative values, these'','//
     :       ''' have been set to zero.'')') NNEG
             CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF
              
         CALL TOTAL(MEMD(PSFPNT(N)),NX,NY,PT)
         IF(DABS(PT-1.0D0).GT.1.0E-5) THEN
            WRITE(CHARS,
     :       '(''! Warning, PSF '',I4,'//
     :    ''' totals '',F15.7,'' - renormalising.'')') N,PT
            CALL UMSPUT(CHARS,1,0,ISTAT)
            CALL MULC(MEMD(PSFPNT(N)),NX,NY,1.0D0/PT,
     :                MEMD(PSFPNT(N)))
         ENDIF
              
         CALL DFILL(MEMD(PSFPNT(N)),NX,NY,MEMD(PSFFFT(N)))
         CALL DFOURT(MEMD(PSFFFT(N)),DIMS,2,-1,0,SCRATCH)
                                                         
C At this point we have finished with the PSF so it can be closed
C and the memory freed - this is NOT the case with modifyable
C PSFs where we will need to change them later
         CALL UIMCLO(IDP(N),ISTAT)
         CALL UDMFRE(PSFPNT(N),7,ISTAT)
      ENDDO
           
C Create other arrays needed, workspace etc
      CALL UDMGET(NX*NY*2,7,WORK,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :      '! Unable to allocate memory for workspace array',
     :               1,0,ISTAT)
         GO TO 99
      ENDIF
             
C REST is the total restored array - stars+continuum
      CALL UDMGET(NX*NY,7,REST,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :      '! Unable to allocate memory for restored array',
     :               1,0,ISTAT)
         GO TO 99
      ENDIF
             
C PS is the continuum part
      CALL UDMGET(NX*NY,7,PS,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :      '! Unable to allocate memory for continuum array',
     :               1,0,ISTAT)
         GO TO 99
      ENDIF

C PSM is the point sources array(s) (only if there are any)
      IF(NSTARS.GT.0) THEN
         DO N=1,NPSF
            CALL UDMGET(NX*NY,7,PSM(N),ISTAT)
            IF(ISTAT.NE.0) THEN
               CALL UMSPUT(
     :         '! Unable to allocate memory for PSM array',
     :                  1,0,ISTAT)
               GO TO 99
            ENDIF
         ENDDO
      ENDIF

C CH is the 'entropy' array
      CALL UDMGET(NX*NY,7,CH,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :      '! Unable to allocate memory for CH array',
     :               1,0,ISTAT)
         GO TO 99
      ENDIF

C PHS 
      CALL UDMGET(NX*NY,7,PHS,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT(
     :      '! Unable to allocate memory for PHS array',
     :              1,0,ISTAT)
        GO TO 99
      ENDIF
              
C Get space for and initialise the correction factor array
      CALL UDMGET(NX*NY,7,CF,ISTAT)
      IF(ISTAT.NE.0) THEN
         CALL UMSPUT('! Unable to create corr. fac. image',
     :                    1,0,ISTAT)
         GO TO 99
      ENDIF
            
C Prepare the first estimate - just a flat image with 
C the fluxes of the point sources deduced from the third column
C of the input list, using the magzero parameter to scale them.
C The remaining flux goes into the background.
C 
C If there are no specified star positions
C we can just fill the array with a constant
      IF(NSTARS.GT.0) THEN

C First work out the total of the star intensities
         ST=0.0D0
         DO I=1,NSTARS
            IF(STARS(4,I).NE.0) 
     :      ST=ST+STARS(3,I)
         ENDDO

C Write out the initial background value and check that it is
C positive
         IF(VERBOSE) THEN
          WRITE(CHARS,
     :   '(''-Setting background in first estimate to: '','//
     :                 'F12.3)') (T-ST)/DFLOAT(NX*NY)
          CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

C Check that the background would be positive - if not set it
C to zero and get it renormalised later
         IF((T-ST)/DFLOAT(NX*NY).LT.0.0D0) THEN
            CALL UMSPUT(
     :     '! Error - negative initial background, setting to zero',
     :                  1,0,ISTAT)
            CALL FILCON(MEMD(PS),NX,NY,0.0D0)
         ELSE
            CALL FILCON(MEMD(PS),NX,NY,(1.0D0-ST/T)/DFLOAT(NX*NY))
         ENDIF

         DO N=1,NPSF
            CALL FILCON(MEMD(PSM(N)),NX,NY,0.0D0)

            CALL MKINIM(MEMD(PSM(N)),
     :               N,NX,NY,STARS,NSTARS,
     :               1.0D0/T)
         ENDDO
      ELSE
         CALL FILCON(MEMD(PS),NX,NY,1.0D0/DFLOAT(NX*NY))
      ENDIF
           
C Now we can start the iterative procedure
      CALL UCLGSI('niter',NITER,ISTAT)
                                      
      DO N=1,NITER
                  
         IF(VERBOSE) THEN
          CALL UMSPUT(' ',1,0,ISTAT)
            WRITE(CHARS,'(''# Starting iteration '',I5,'//
     :            '''/'',I5,''-----------'')') N,NITER
          CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

C Check for non-negativity
         CALL ZAPNEG(MEMD(PS),NX,NY,NNEG,VMIN,0.0D0)
         IF(NNEG.GT.0) THEN
            WRITE(CHARS,'(''! Background image '','//
     :         '''contains '',I7,'' negative values'')') NNEG
            CALL UMSPUT(CHARS,1,0,ISTAT)
            WRITE(CHARS,'(''--Minimum value: '',D20.10)') VMIN
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

         IF(NSTARS.GT.0) THEN
          DO I=1,NPSF
            CALL ZAPNEG(MEMD(PSM(I)),NX,NY,NNEG,VMIN,0.0D0)
            IF(NNEG.GT.0) THEN
               WRITE(CHARS,'(''! Point Source image '',I3,'//
     :         ''' contains '',I7,'' negative values'')') I,NNEG
               CALL UMSPUT(CHARS,1,0,ISTAT)
               WRITE(CHARS,'(''--Minimum value: '',D20.10)') VMIN
               CALL UMSPUT(CHARS,1,0,ISTAT)
            ENDIF
          ENDDO
         ENDIF
                                                 
C Check that everything is correctly normalised
         CALL TOTAL(MEMD(PS),NX,NY,TT)

         IF(NSTARS.GT.0) THEN
            ST=0.0D0
            DO I=1,NPSF
               CALL TOTAL(MEMD(PSM(I)),NX,NY,TEM)
               ST=ST+TEM
               TT=TT+TEM
            ENDDO

            WRITE(CHARS,'(''--Flux distribution, Stars: '','//
     :    'F6.2,''% Background: '',F6.2,''%.'')') 100.0D0*ST,
     :    100.0D0*(1.0D0-ST)
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

         IF(DABS(TT-1.0D0).GT.1.0D-10) THEN
            WRITE(CHARS,'('//
     :     '''! Warning, sum of predicted image is '',D20.10,'//
     :     '''(not renormalised).'')') TT
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

C First add the continuum and point source arrays (if the latter
C exists). Only the first (centered) point source array.
         IF(NSTARS.GT.0) THEN
            CALL ADDIM(MEMD(PS),MEMD(PSM(1)),NX,NY,MEMD(REST))

C Convolve the current estimate with the PSF to give PHS     
            CALL DCONV(MEMD(REST),NX,NY,MEMD(WORK),
     :                 MEMD(PSFFFT(1)),MEMD(PHS),1)

C Now convolve the rest of the PSFs and add those in too
            IF(NPSF.GT.1) THEN
               DO I=2,NPSF
                  CALL DCONV(MEMD(PSM(I)),NX,NY,MEMD(WORK),
     :                 MEMD(PSFFFT(I)),MEMD(W),1)
                  CALL ADDIM(MEMD(W),MEMD(PHS),
     :                 NX,NY,MEMD(PHS))
               ENDDO
            ENDIF
         ELSE
           CALL DCONV(MEMD(PS),NX,NY,MEMD(WORK),
     :                 MEMD(PSFFFT(1)),MEMD(PHS),1)
         ENDIF
              
C If we have a floating prior we now convolve with the smoothing
C kernel to get the default solution
C
C NB, this is the version 0.5D bugfix - the ,1 at the
C end was missing in earlier versions.
         IF(FPRIOR) THEN
            CALL DCONV(MEMD(PS),NX,NY,MEMD(WORK),
     :                 MEMD(FPFFT),MEMD(PHP),1)
         ENDIF

C Calculate and display the residual information to assist closeness
C of fit assessment
         IF(VERBOSE) THEN
             CALL RESINF(MEMD(PHS),MEMD(DATPNT),
     :                    NX,NY,RMSRES,RESMAX,IRMAX,JRMAX)
               
          WRITE(CHARS,'(''--RMS residual: '',D10.4,'//
     :                 ''' Max residual of '',D10.4,'//
     :                 ''' at ('',I4,'','',I4,'')'')')
     :             RMSRES*T,RESMAX*T,IRMAX,JRMAX
          CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF
              
C Calculate the sum of the continuum component
         CALL TOTAL(MEMD(PS),NX,NY,SMB)

C If we are sub-sampling we must do so here.
         CALL TOTAL(MEMD(PHS),NX,NY,SMB)

         IF(SUBSAM) CALL REBIN(MEMD(PHS),NX,NY,NBX,NBY)

         CALL TOTAL(MEMD(PHS),NX,NY,SMB)

C Divide the data by this convolution
         CALL DIVIDE(MEMD(DATPNT),MEMD(PHS),NX,NY,MEMD(W))

C If there is a floating prior we need to calculate the RP
C array
         IF(FPRIOR) THEN
            CALL DIVIDE(MEMD(PS),MEMD(PHP),NX,NY,MEMD(RP))
              
C Calculate the entropy and likelihood and fill the array CH
C This is a bit different if there is a floating prior so we
C use a different routine
            CALL FNTROP(MEMD(PS),MEMD(PHS),MEMD(DATPNT),MEMD(PHP),
     :               NX,NY,MEMD(CH),SS,HH)
         ELSE
            CALL ENTROP(MEMD(PS),MEMD(PHS),MEMD(DATPNT),SMB,
     :               NX,NY,MEMD(CH),SS,HH)
         ENDIF

C Calculate the objective function
         QQ=HH+CS*SS

C Some diagnostics
         IF(VERBOSE) THEN
          WRITE(CHARS,'(''--Log.Lik: '','//
     :         'E15.6,'' Ent: '',E15.6,'//
     :         ''' Obj.Func: '',E20.12)')
     :    HH,SS,QQ
          CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

         IF(QQ.LT.LV) THEN
            WRITE(CHARS,
     :  '(''! Warning, objective function fell by '',D15.7)')
     :      LV-QQ
            CALL UMSPUT(CHARS,1,0,ISTAT)

            IF(.NOT.ACCEL) THEN
               CALL UMSPUT('! Turning on acceleration.',
     :            1,0,ISTAT)

C To avoid catastrophic instabilities we now need to turn on 
C acceleration. First we need to allocate the space we will need
               ACCEL=.TRUE.

C This one is only needed for the floating prior case
               IF(FPRIOR) THEN
                  CALL UDMGET(NX*NY,7,DPHP,ISTAT)
                  IF(ISTAT.NE.0) THEN
                     CALL UMSPUT(
     :         '! Unable to allocate memory for DPHP array',
     :                       1,0,ISTAT)
                     GO TO 99
                  ENDIF
               ENDIF

               CALL UDMGET(NX*NY,7,DPHS,ISTAT)
               IF(ISTAT.NE.0) THEN
                  CALL UMSPUT(
     :           '! Unable to allocate memory for DPHS array',
     :                    1,0,ISTAT)
                  GO TO 99
               ENDIF
            ENDIF
         ELSE
            IF(VERBOSE.AND.N.GT.1) THEN
               WRITE(CHARS,
     :  '(''--Objective function increased by '',D15.7)')
     :         QQ-LV
               CALL UMSPUT(CHARS,1,0,ISTAT)
            ENDIF
         ENDIF
         LV=QQ

C Convolve again, this time by the rotated PSF to get the correction
C factor array
         CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :          MEMD(PSFFFT(1)),MEMD(CF),-1)
                                              
C If we have a floating prior we need to calculate another
C correction factor array (this is also a correlation)
         IF(FPRIOR) THEN
            CALL DCONV(MEMD(RP),NX,NY,MEMD(WORK),
     :          MEMD(FPFFT),MEMD(CFP),-1)
         ENDIF

C If we have acceleration enabled we should now calculate a couple
C more arrays, even more convolutions unfortunately...
         IF(ACCEL) THEN
               
C First calculate DPSS (into MEMD(W)):
C (different code for floating or stationary prior)
             IF(FPRIOR) THEN
                CALL FIDPSS(MEMD(PS),MEMD(PSM(1)),
     :               MEMD(CF),MEMD(CFP),
     :            MEMD(CH),NX,NY,CS,SS,NSTARS,MEMD(W))
             ELSE
                CALL CIDPSS(MEMD(PS),MEMD(PSM(1)),
     :          MEMD(CF),MEMD(CH),NX,NY,CS,SMB,SS,NSTARS,
     :          MEMD(W))
             ENDIF

C Convolve to give DPHS
             CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :           MEMD(PSFFFT(1)),MEMD(DPHS),1)

C Maybe we should rebin here?
             IF(SUBSAM) CALL REBIN(MEMD(DPHS),NX,NY,NBX,NBY)

C Calculate DPS (into MEMD(W) again):
C If there are no stars this can be ignored as W
C already has the correct contents
             IF(NSTARS.GT.0) 
     :           CALL FILDPS(MEMD(W),MEMD(PSM(1)),MEMD(CF),
     :           NX,NY,CS,SS,MEMD(W))

C Convolve this to give DPHP (floating prior only):
             IF(FPRIOR) THEN
                CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :               MEMD(FPFFT),MEMD(DPHP),1)
               
C Calculate the acceleration factor
                CALL FACCL(MEMD(DATPNT),MEMD(PS),MEMD(W),
     :           MEMD(PHS),MEMD(DPHS),MEMD(PHP),MEMD(DPHP),
     :           NX,NY,CS,XL,VERBOSE)
             ELSE
                CALL ACCL(MEMD(DATPNT),MEMD(PS),MEMD(W),
     :           MEMD(PHS),MEMD(DPHS),
     :           NX,NY,CS,XL,VERBOSE)
             ENDIF

             IF(VERBOSE) THEN
              WRITE(CHARS,
     ;       '(''--Using acceleration factor of: '',F10.4)') XL
              CALL UMSPUT(CHARS,1,0,ISTAT)
             ENDIF
         ENDIF

C Apply the correction, first to the continuous image
         IF(FPRIOR) THEN
            CALL FORCON(MEMD(PS),MEMD(CF),MEMD(CFP),MEMD(PSM(1)),
     :                  MEMD(CH),NX,NY,XL,CS,SS,TST)
         ELSE
            CALL CORCON(MEMD(PS),MEMD(CF),
     :                  MEMD(CH),NX,NY,XL,CS,SS,SMB,TST)
         ENDIF

         IF(VERBOSE) THEN
            WRITE(CHARS,
     :       '(''--Convergence parameter: '',D20.8)') TST
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF

C and to the point source image (if there is one)
         IF(NSTARS.GT.0) THEN  
            CALL CORSTR(MEMD(PSM(1)),MEMD(CF),1,NX,NY,
     :             XL,T,STARS,NSTARS,SS,CS,MAGZERO,VERBOSE)

C now if there are more PSFs we need to do this again
            IF(NPSF.GT.1) THEN
              DO I=2,NPSF
               CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :             MEMD(PSFFFT(I)),MEMD(CF),-1)
               CALL CORSTR(MEMD(PSM(I)),MEMD(CF),I,NX,NY,
     :            XL,T,STARS,NSTARS,SS,CS,MAGZERO,VERBOSE)
              ENDDO
            ENDIF
         ENDIF
      ENDDO
           
C Add up the components of the result (if there are 
C any point sources)
      IF(NSTARS.GT.0) THEN
         CALL ADDIM(MEMD(PS),MEMD(PSM(1)),NX,NY,MEMD(REST))

         IF(NPSF.GT.1) THEN
            DO N=2,NPSF
               CALL ADDIM(MEMD(REST),MEMD(PSM(N)),NX,NY,MEMD(REST))
            ENDDO
         ENDIF

C Multiply up by the total in the input array (un-normalise)
         CALL MULC(MEMD(REST),NX,NY,T,MEMD(REST))
      ELSE
         CALL MULC(MEMD(PS),NX,NY,T,MEMD(REST))
      ENDIF

C Write out the results: 
C First the deconvolved image
      IF(VERBOSE)
     :  CALL UMSPUT('-Writing output deconvolved image: '//DECON
     :              ,1,0,ISTAT)
                               
      CALL UIPS2D(IDDC,1,NX,1,NY,MEMD(REST),ISTAT)
      CALL UIMCLO(IDDC,ISTAT)
                              
C Then the background image if requested
      IF(BACK.NE.' ') THEN
         IF(VERBOSE)
     :     CALL UMSPUT(
     :       '-Writing output smooth background image: '//BACK
     :                 ,1,0,ISTAT)
                               
         CALL MULC(MEMD(PS),NX,NY,T,MEMD(PS))
         CALL UIPS2D(IDPS,1,NX,1,NY,MEMD(PS),ISTAT)
         CALL UIMCLO(IDPS,ISTAT)
      ENDIF

C And now the text list
      IF(OUTTAB.NE.' '.AND.NSTARS.GT.0) THEN

C First we write in some information which may be of interest - mainly
C to record how the file was produced
         WRITE(8,'(''# Created by PLUCY V0.70'')')
         WRITE(8,'(''#'')')
         WRITE(8,'(''# Input data image: '',A30)') DATA
         WRITE(8,'(''# Input PSF image(s): '',A50)') PSFS
         WRITE(8,'(''# Fractional entropy strength: '',F10.4)') CS
         WRITE(8,'(''# Number of interations: '',I6)') NITER
         WRITE(8,'(''# Output deconvolved image: '',A50)') DECON
         WRITE(8,'(''# Input star list: '',A50)') FILE
         WRITE(8,'(''# Floating prior was used.'')')
         WRITE(8,'(''# Sigma of smoothing kernel: '',F10.3)') GW

         IF(ACCEL) THEN
            WRITE(8,'(''# Acceleration was used.'')')
         ELSE
            WRITE(8,'(''# Acceleration was not used.'')')
         ENDIF

         WRITE(8,'(''# Magnitude zero point: '',F10.3)') MAGZERO

         WRITE(8,
     :   '(''# Number  X        Y       InputFlux'','//
     :   '''   MeasuredFlux MeasuredMag'')')

         CALL WRITAB(STARS,NSTARS,MAGZERO,NPSF)
         IF(VERBOSE)
     :     CALL UMSPUT(
     :       '-Writing output photometry list: '//OUTTAB
     :                 ,1,0,ISTAT)
      ENDIF

99    CONTINUE
               
C Close the text file
      CLOSE(8)

C Close all the images
      CALL UIMCLO(IDD,ISTAT)
           
C Free all the dynamic arrays
      CALL UDMFRE(DATPNT,7,ISTAT)

      DO I=1,NPSF
         CALL UDMFRE(PSFFFT(I),7,ISTAT)
      ENDDO

      CALL UDMFRE(PHS,7,ISTAT)
      CALL UDMFRE(PS,7,ISTAT)

      IF(NSTARS.GT.0) THEN
         DO I=1,NPSF
            CALL UDMFRE(PSM(I),7,ISTAT)
         ENDDO
      ENDIF

      CALL UDMFRE(CF,7,ISTAT)
      CALL UDMFRE(WORK,7,ISTAT)
      CALL UDMFRE(REST,7,ISTAT)
      CALL UDMFRE(W,7,ISTAT)

C Another one if acceleration is on
      CALL UDMFRE(DPHS,7,ISTAT)
      
C Extra ones if there is a floating prior
      IF(FPRIOR) THEN
        CALL UDMFRE(PHP,7,ISTAT)
        CALL UDMFRE(FPFFT,7,ISTAT)
        CALL UDMFRE(RP,7,ISTAT)
        CALL UDMFRE(CFP,7,ISTAT)

C Even more if there is acceleration also
        IF(ACCEL) THEN
           CALL UDMFRE(DPHP,7,ISTAT)
        ENDIF
      ENDIF
                                      
C Close down MIDAS
CM    CALL STSEPI    !M
                       
      END
         
      SUBROUTINE MULC(IN1,NX,NY,F,OUT)
C
C Just multiply one array by a constant
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J
      DOUBLE PRECISION IN1(NX,NY),OUT(NX,NY),F
                                              
      DO J=1,NY
         DO I=1,NX
            OUT(I,J)=IN1(I,J)*F
         ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE DIVIDE(IN1,IN2,NX,NY,OUT)
C
C Just divide one array by another of the same size
C If the denominator is zero the result is set to zero also.
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J
      DOUBLE PRECISION IN1(NX,NY),IN2(NX,NY),OUT(NX,NY)
      INTEGER NZERO,ISTAT
      CHARACTER*80 CHARS

      NZERO=0
                                                       
      DO J=1,NY
         DO I=1,NX
            IF(IN2(I,J).EQ.0.0D0) THEN
               OUT(I,J)=0.0D0
               NZERO=NZERO+1
            ELSE
               OUT(I,J)=IN1(I,J)/IN2(I,J)
            ENDIF
         ENDDO
      ENDDO
           
      IF(NZERO.GT.0) THEN
         WRITE(CHARS,
     : '(''! '',i10,'' attempts to divide by zero.'')') NZERO
         CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF

      RETURN
      END
         
      SUBROUTINE TOTAL(DATA,NX,NY,TOT)
C
C Add up all elements in an image
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J
      DOUBLE PRECISION DATA(NX,NY),TOT
                                      
      TOT=0.0D0
      DO J=1,NY
         DO I=1,NX
             TOT=TOT+DATA(I,J)
         ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE MKINIM(IMAGE,NPSF,NX,NY,STARS,NSTARS,T)
C
C Make a suitable starting image from the list of Stars
C provided. The values are those given in the 3rd column
C of the STARS list multiplied by the value of T.
C
C It is assumed that the background is empty on entry.
C

      IMPLICIT NONE

      INTEGER NX,NY,NSTARS
      DOUBLE PRECISION IMAGE(NX,NY),T,STARS(5,NSTARS)
      INTEGER I,NPSF,ISTAT
      CHARACTER*80 CHARS

      DO I=1,NSTARS
         IF(STARS(3,I).NE.0.0D0) THEN
           IF(NINT(STARS(4,I)).EQ.NPSF) THEN
              IF(IMAGE(NINT(STARS(1,I)),NINT(STARS(2,I)))
     :           .NE.0.0D0) THEN
                 WRITE(CHARS,
     : '(''! Pixel ('',I4,'','',I4,'') already has a star in it.'')')
     :          NINT(STARS(1,I)),NINT(STARS(2,I))
                CALL UMSPUT(CHARS,1,0,ISTAT)
                STARS(4,I)=0.0D0
              ELSE 
                 IMAGE(NINT(STARS(1,I)),NINT(STARS(2,I)))=
     :          T*STARS(3,I)
              ENDIF
           ENDIF
         ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE CORCON(PS,CF,CH,NX,NY,XL,CS,SS,SMB,TST)
C
C Apply the correction to the continuous background image
C including an entropy term. Constant prior version.
C
C This routine also calculates values for a normalisation
C and convergence test.
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION PS(NX,NY),CF(NX,NY)
      DOUBLE PRECISION CH(NX,NY),CS,SS,SMB,XL
      DOUBLE PRECISION DH,DS,DT,TST,VH,VS,XMAX
      INTEGER I,J

      DH=0.0D0
      DS=0.0D0
      XMAX=1.0D30
      DT=0.0D0

      DO J=1,NY
         DO I=1,NX
            VH=PS(I,J)*(CF(I,J)-1.0D0)
            VS=CS*(-PS(I,J)*SS/SMB-CH(I,J))
            DH=DH+VH*VH
            DS=DS+VS*VS
            DT=DT+(VH+VS)**2
            PS(I,J)=PS(I,J)+XL*(VH+VS)
         ENDDO
      ENDDO

      TST=DSQRT(DT)/(DSQRT(DH)+DSQRT(DS))

      RETURN
      END

      SUBROUTINE FORCON(PS,CF,CFP,PSM,
     :           CH,NX,NY,XL,CS,SS,TST)
C
C Apply the correction to the continuous background image
C including an entropy term. Floating prior version.
C
C This routine also calculates values for a normalisation
C and convergence test.
C
C Modified for the new, simpler normalisation scheme and also
C removed the SMB argument - Richard Hook, Dec 14th 1994
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION PS(NX,NY),CF(NX,NY),CFP(NX,NY)
      DOUBLE PRECISION CH(NX,NY),PSM(NX,NY),CS,SS,XL
      DOUBLE PRECISION DH,DS,DT,TST,VH,VS,SSS
      INTEGER I,J

      DH=0.0D0
      DS=0.0D0
      DT=0.0D0

      SSS=0.0D0

      DO J=1,NY
         DO I=1,NX
            VH=PS(I,J)*(CF(I,J)-1.0D0)
            VS=CS*PS(I,J)*(CFP(I,J)-1.0D0)+
     :         CS*(-PS(I,J)*SS-CH(I,J))

            SSS=SSS+VS+VH
            DH=DH+VH*VH
            DS=DS+VS*VS
            DT=DT+(VH+VS)**2
            PS(I,J)=PS(I,J)+XL*(VH+VS)
         ENDDO
      ENDDO

      TST=DSQRT(DT)/(DSQRT(DH)+DSQRT(DS))

      RETURN
      END

      SUBROUTINE CORSTR(PSM,CF,NPSF,NX,NY,XL,T,STARS,NSTARS,
     :                  SS,CS,MAGZERO,VERBOSE)
C
C Apply the correction to the point sources - they have no
C entropy term and are just multiplied by the correction
C factor
C
C Modified to include normalisation changes, Dec 94
C
      IMPLICIT NONE

      INTEGER NX,NY,NSTARS,ISTAT,NPSF
      DOUBLE PRECISION PSM(NX,NY),CF(NX,NY),SS,CS
      DOUBLE PRECISION STARS(5,NSTARS),T,XL,MAGZERO
      CHARACTER*80 CHARS
      LOGICAL VERBOSE
     
      INTEGER I,J,IX,IY

      DO J=1,NY
         DO I=1,NX
            PSM(I,J)=PSM(I,J)*(1.0D0+XL*(CF(I,J)-1.0D0-CS*SS))
         ENDDO
      ENDDO

      DO I=1,NSTARS
        IF(STARS(3,I).NE.0.0D0 .AND.
     :     NINT(STARS(4,I)).EQ.NPSF) THEN
         IX=NINT(STARS(1,I))
         IY=NINT(STARS(2,I))
         STARS(5,I)=PSM(IX,IY)*T

         IF(I.LT.5.AND.VERBOSE.AND.STARS(5,I).GT.0.0D0) THEN
          IF(NPSF.GT.1) THEN
           WRITE(CHARS,'(''--Flux in *'',I4,'' is: '','//
     :         'G14.6,'' Mag: '',F9.3,'' (PSF: '',I4,'')'')') 
     :         I,PSM(IX,IY)*T,
     :         MAGZERO-2.5*DLOG10(PSM(IX,IY)*T),NPSF
           ELSE
           WRITE(CHARS,'(''--Flux in *'',I4,'' is: '','//
     :         'G14.6,'' Mag: '',F9.3)') I,PSM(IX,IY)*T,
     :         MAGZERO-2.5*DLOG10(PSM(IX,IY)*T)
           ENDIF
           CALL UMSPUT(CHARS,1,0,ISTAT)
          ENDIF
         ENDIF
         IF(I.EQ.6.AND.VERBOSE) 
     :        CALL UMSPUT('    Etc,Etc...',1,0,ISTAT)
      ENDDO

      RETURN
      END

      SUBROUTINE WRITAB(STARS,NSTARS,MAGZERO,NPSF)
C
C Write out the final results as an ASCII table
C It is assumed that this file has already been opened
C on FORTRAN unit 8.
C
      IMPLICIT NONE

      INTEGER NSTARS
      DOUBLE PRECISION STARS(5,NSTARS),MAGZERO

      INTEGER I,IX,IY,K,NPSF

      DO I=1,NSTARS
         IX=NINT(STARS(1,I))
         IY=NINT(STARS(2,I))

         IF(STARS(3,I).NE.0.0D0 .AND. STARS(5,I).GT.0.0D0) THEN
          IF(NPSF.GT.1) THEN
           WRITE(8,'(I6,1X,F7.2,1X,F7.2,1X,'//
     :         'F12.2,1X,F12.2,1X,F7.3,1X,I5)')
     :         I,(STARS(K,I),K=1,3),STARS(5,I),
     :         MAGZERO-2.5*DLOG10(STARS(5,I)),
     :         NINT(STARS(4,I))
          ELSE
           WRITE(8,'(I6,1X,F7.2,1X,F7.2,1X,'//
     :         'F12.2,1X,F12.2,1X,F7.3)')
     :         I,(STARS(K,I),K=1,3),STARS(5,I),
     :         MAGZERO-2.5*DLOG10(STARS(5,I))
          ENDIF
         ELSE
          IF(NPSF.GT.1) THEN
           WRITE(8,'(I6,1X,F7.2,1X,F7.2,1X,'//
     :         'F12.2,1X,F12.2,1X,F7.3,1X,I5)')
     :         I,(STARS(K,I),K=1,3),0.0D0,
     :         -100.0,NINT(STARS(4,I))
          ELSE
           WRITE(8,'(I6,1X,F7.2,1X,F7.2,1X,'//
     :         'F12.2,1X,F12.2,1X,F7.3)')
     :         I,(STARS(K,I),K=1,3),0.0D0,-100.0
          ENDIF
         ENDIF

      ENDDO

      END

      SUBROUTINE MULIM(A,B,NX,NY,OUT)
C
C Multiply an image by another
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      DOUBLE PRECISION A(NX,NY),B(NX,NY),OUT(NX,NY)

      DO J=1,NY
         DO I=1,NX
            OUT(I,J)=A(I,J)*B(I,J)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE ADDIM(A,B,NX,NY,OUT)
C
C Add an image to another
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      DOUBLE PRECISION A(NX,NY),B(NX,NY),OUT(NX,NY)

      DO J=1,NY
         DO I=1,NX
            OUT(I,J)=A(I,J)+B(I,J)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE ENTROP(PS,PHS,PHT,SMB,NX,NY,CH,SS,HH)
C
C Calculate the entropy and likelihood
C as well as filling the array CH. This
C version is for the stationary prior.
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION PS(NX,NY),PHS(NX,NY),PHT(NX,NY)
      DOUBLE PRECISION SMB,CH(NX,NY),SS,HH

      INTEGER I,J

      SS=0.0D0
      HH=0.0D0

      DO J=1,NY
         DO I=1,NX
            CH(I,J)=0.0D0
            IF(PS(I,J).GT.0.0D0)
     :        CH(I,J)=PS(I,J)/SMB*DLOG(PS(I,J)/SMB)
            SS=SS-CH(I,J)
            IF(PHS(I,J).GT.0.0D0) HH=HH+PHT(I,J)*
     :                              DLOG(PHS(I,J))
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE FNTROP(PS,PHS,PHT,PHP,NX,NY,CH,SS,HH)
C
C Calculate the entropy and likelihood
C as well as filling the array CH
C
C This version uses a floating prior
C
C NB - modified in version 0.6 for a simpler normalisation
C        Richard Hook, August 1994
C
C Note also that I use CH and Leon uses SF:
C
C CH(i,j) = PS(i,j)*SF(i,j)
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION PS(NX,NY),PHS(NX,NY),PHT(NX,NY)
      DOUBLE PRECISION CH(NX,NY),PHP(NX,NY),SS,HH

      INTEGER I,J

      SS=0.0D0
      HH=0.0D0

      DO J=1,NY
         DO I=1,NX
            CH(I,J)=0.0D0
            IF(PS(I,J).GT.0.0D0 .AND. PHP(I,J).GT.0.0D0)
     :        CH(I,J)=PS(I,J)*DLOG(PS(I,J)/PHP(I,J))
            SS=SS-CH(I,J)
            IF(PHS(I,J).GT.0.0D0) HH=HH+PHT(I,J)*
     :                              DLOG(PHS(I,J))
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE PNTROP(RHO,RHOP,NX,NY,TF,TT)
C
C Calculate the entropy for the PSF in the case
C of a regularized, modifyable PSF
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION RHO(NX,NY),RHOP(NX,NY)
      DOUBLE PRECISION TF(NX,NY),TT

      INTEGER I,J

      TT=0.0D0

      DO J=1,NY
         DO I=1,NX
            TF(I,J)=0.0D0
            IF(RHO(I,J).GT.0.0D0 .AND. RHOP(I,J).GT.0.0D0)
     :        TF(I,J)=DLOG(RHO(I,J)/RHOP(I,J))
            TT=TT-RHO(I,J)*TF(I,J)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE RESINF(IM1,IM2,NX,NY,RMSRES,RESMAX,
     :                   IRMAX,JRMAX)
C
C Calculate the RMS residual between two images and the
C largest residual.
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,IRMAX,JRMAX
      DOUBLE PRECISION RMSRES,RESMAX,IM1(NX,NY),IM2(NX,NY)
                                                          
      INTEGER I,J
      DOUBLE PRECISION RES,T
                            
      T=0.0D0
      RESMAX=0.0D0
                  
      DO J=1,NY
         DO I=1,NX
            RES=IM1(I,J)-IM2(I,J)
            T=T+RES*RES
            IF(DABS(RES).GT.DABS(RESMAX)) THEN
               RESMAX=RES
               IRMAX=I
               JRMAX=J
            ENDIF
         ENDDO
      ENDDO
           
      RMSRES=DSQRT(T/DBLE(NX*NY))
                                 
      RETURN
      END
         
      SUBROUTINE REBIN(IN,NX,NY,NBX,NBY)
C
C Re-bin an array by replacing groups of pixels by their
C average. This is needed in the sub-sampling case.
C Flux is conserved and the operation is done 'in place'
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,NBX,NBY
      DOUBLE PRECISION IN(NX,NY)
                                
      INTEGER I,J,K,L
      DOUBLE PRECISION T,DN
                           
      DN=DBLE(NBX*NBY)
                      
      DO J=1,NY,NBY
         DO I=1,NX,NBX
            T=0D0
            DO L=1,NBY
               DO K=1,NBX
                  T=T+IN(I+K-1,J+L-1)
               ENDDO
            ENDDO
            DO L=1,NBY
               DO K=1,NBX
                  IN(I+K-1,J+L-1)=T/DN
               ENDDO
            ENDDO
         ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE ZAPNEG(DATA,NX,NY,NNEG,VMIN,VZAP)
C
C Find negative points in an array and set them to VZAP
C Also return the number of negative points and the smallest
C value.
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,NNEG
      DOUBLE PRECISION DATA(NX,NY),VMIN,VZAP
      INTEGER I,J
                 
      NNEG=0
      VMIN=1.0D20
            
      DO J=1,NY
         DO I=1,NX
            IF(DATA(I,J).LT.0.0D0) THEN
               NNEG=NNEG+1
               IF(DATA(I,J).LT.VMIN) VMIN=DATA(I,J)
               DATA(I,J)=VZAP
            ENDIF
         ENDDO
      ENDDO
           
      RETURN
      END

      SUBROUTINE GETSTR(FILE,STARS,NX,NY,NPSF,
     :                  W,MAXSTR,MAGZERO,NSTARS,ISTAT)
C
C Read in a list of free format star positions (in pixels)
C
C Each line has four entries:
C  - the X position
C  - the Y position
C  - some sort of relative brightness expressed as a magnitude
C    relative to MAGZERO
C  - the number of the matching PSF (only if NPSF>1)
C
C Comments may be put in a line after a ! or #, blank lines
C are ignored.
C
      IMPLICIT NONE

      INTEGER MAXSTR,NSTARS,ISTAT,NX,NY,I,NPSF
      DOUBLE PRECISION STARS(5,MAXSTR),MAGZERO
      DOUBLE PRECISION W(NX,NY)
      CHARACTER*80 FILE,CHARS,LINE

      INTEGER N,IOS,NN

      ISTAT=0
      NSTARS=0

C First open the file
      OPEN(12,FILE=FILE,STATUS='OLD',IOSTAT=IOS)
      IF(IOS.NE.0) THEN
         CALL UMSPUT('! Unable to open star positions list file',
     :               1,0,ISTAT)
         ISTAT=1
         RETURN
      ENDIF

      N=0

C Read in the values until there are no more
      DO NN=1,MAXSTR+1000
         IF(N.EQ.MAXSTR+1) THEN
            CALL UMSPUT('! Too many stars!',1,0,ISTAT)
            ISTAT=1
            GO TO 99 
         ENDIF

C Read a line from the file
         READ(12,'(A80)',END=99) LINE 

C Remove comments
         DO I=1,80
            IF(LINE(I:I).EQ.'#' .OR.
     :         LINE(I:I).EQ.'!') THEN
               LINE(I:80)=' '
            ENDIF
         ENDDO

         IF(LINE.NE.' ') THEN
            N=N+1

C If there is only one PSF we don't try to get the PSF
C number, we just default it to 1
            IF(NPSF.GT.1) THEN
             READ(LINE,*,IOSTAT=IOS) (STARS(I,N),I=1,4)
            ELSE
             READ(LINE,*,IOSTAT=IOS) (STARS(I,N),I=1,3)
      write(chars,'('' Stars: '',i5,3f12.2)') n,(STARS(I,N),I=1,3)

      call umsput(chars,1,0,istat)
             STARS(4,N)=1.0D0
            ENDIF
 
            IF(IOS.NE.0) THEN
               CALL UMSPUT('! Error reading StarList file',
     :                     1,0,ISTAT)
               GO TO 99
            ENDIF

C Convert from magnitudes to fluxes
            STARS(3,N)=10.0**(-0.4*(STARS(3,N)-MAGZERO))

C Check for a star already there
            IF(W(NINT(STARS(1,N)),NINT(STARS(2,N))).NE.
     :         0.0D0) THEN
        write(chars,'('' W val: '',f12.4)') 
     :     W(NINT(STARS(1,N)),NINT(STARS(2,N)))
        call umsput(chars,1,0,istat)
               STARS(3,N)=0.0D0
               STARS(4,N)=0.0D0
               WRITE(CHARS,'(''! Warning duplicate star at '','//
     :          '''X: '',I4,'' Y: '',I4)') NINT(STARS(1,N)),
     :          NINT(STARS(2,N))
               CALL UMSPUT(CHARS,1,0,ISTAT)
            ENDIF

C Check for out of bounds stars
            IF(NINT(STARS(1,N)).LT.1 .OR.
     :         NINT(STARS(1,N)).GT.NX .OR.
     :         NINT(STARS(2,N)).LT.1 .OR.
     :         NINT(STARS(2,N)).GT.NY) THEN
               STARS(3,N)=0.0D0
               STARS(4,N)=0.0D0
               WRITE(CHARS,'(''! Warning, star '',I5,'//
     :           ''' is outside image, ignored.'')') N
               CALL UMSPUT(CHARS,1,0,ISTAT)
            ELSE
               W(NINT(STARS(1,N)),NINT(STARS(2,N)))=1.0D0
            ENDIF

            IF(NINT(STARS(4,N)).LT.0 .OR.
     :         NINT(STARS(4,N)).GT.NPSF) THEN
               STARS(3,N)=0.0D0
               STARS(4,N)=0.0D0
               WRITE(CHARS,'(''! Warning, star '',I5,'//
     :           ''' has invalid PSF number, ignored.'')') N
               CALL UMSPUT(CHARS,1,0,ISTAT)
            ENDIF

         ENDIF

         NSTARS=N
      ENDDO

99    CONTINUE
      CLOSE(12)

      RETURN
      END

      SUBROUTINE FGAUSS(IMAGE,NX,NY,GW)
C
C Fill the centre of an empty image with a 2d circular
C gaussian of sigma=GW. We go out as far as 5 sigma.
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION IMAGE(NX,NY),GW

C Local Variables etc
      DOUBLE PRECISION R2,S5S,T,DY,DX
      INTEGER I,J

      S5S=25.0D0*GW*GW
      T=0.0D0

      DO J=1,NY
         DY=DBLE(J-NY/2)
         DO I=1,NX
            DX=DBLE(I-NX/2)
            R2=DX*DX+DY*DY
            IF(R2.LT.S5S) THEN
               IMAGE(I,J)=EXP(-R2/(2.0D0*GW*GW))
               T=T+IMAGE(I,J)
            ELSE
               IMAGE(I,J)=0.0D0
            ENDIF
         ENDDO
      ENDDO

      DO J=1,NY
         DO I=1,NX
            IMAGE(I,J)=IMAGE(I,J)/T
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE FACCL(PHT,PS,DPS,PHS,DPHS,PHP,DPHP,
     :                NX,NY,CS,XL,VERBOSE)
C
C Calculate an optimum acceleration factor for point+smooth
C background restorations. This version has a floating prior.
C
C Based on a routine by Leon Lucy.
C
C Modified to use the simpler normalisation scheme, Dec 1994
C
C Removed the MODPSF flag, Mar 1998
C
      IMPLICIT NONE

      INTEGER NX,NY

      DOUBLE PRECISION PS(NX,NY),DPS(NX,NY)
      DOUBLE PRECISION PHS(NX,NY),DPHS(NX,NY)
      DOUBLE PRECISION PHP(NX,NY),DPHP(NX,NY)
      DOUBLE PRECISION PHT(NX,NY)
      DOUBLE PRECISION CS,XL,XLM1,XLM2,XLM3,SS,HH
      DOUBLE PRECISION XLP,XLU,DQ1,DQ2,QQ,XLM,PSD,PHSD
      DOUBLE PRECISION SS1,SS2,BB,PHPD,ET,TS,YL
      LOGICAL VERBOSE
      INTEGER ISTAT,I,J,NIT
      CHARACTER*80 CHARS

C First calculate the maximum possible acceleration within the 
C positivity restraint
      CALL MAXF(PS,DPS,NX,NY,XLM1)

C and now the same thing in the image plane
      CALL MAXF(PHS,DPHS,NX,NY,XLM2)

C and also for the floating prior
      CALL MAXF(PHP,DPHP,NX,NY,XLM3)

      XL=1.0D0

C We now need a 'safe' value for the acceleration - less than 
C the max permitted but still positive...
      XLM=DMIN1(XLM1,XLM2,XLM3)

      IF(VERBOSE) THEN
       WRITE(CHARS,
     : '(''--Max. acc. possible without introducing neg. values: '','//
     :      'F12.4)') XLM
       CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF

      IF(XLM.GT.1.0D0) THEN
         XLU=1.0D0+0.9D0*(XLM-1.0D0)
      ELSE
         XLU=0.9*XLM
      ENDIF

      NIT=0

C Return here when iterating
45    CONTINUE

      NIT=NIT+1
      XLP=XL

C Derive the first and second derivatives w.r.t. the
C acceleration factor
      SS=0.0D0
      SS1=0.0D0
      SS2=0.0D0
      HH=0.0D0

      DO J=1,NY
         DO I=1,NX
           PSD=PS(I,J)+XL*DPS(I,J)
           PHSD=PHS(I,J)+XL*DPHS(I,J)
           PHPD=PHP(I,J)+XL*DPHP(I,J)

           IF(PHPD.GT.0.0D0 .AND. PSD.GT.0.0D0
     :         .AND. PHSD.GT.0.0D0) THEN
            YL=DLOG(PSD/PHPD)
            ET=DPHP(I,J)/PHPD
            TS=DPS(I,J)/PSD

            SS=SS-PSD*YL
            HH=HH+PHT(I,J)*DLOG(PHSD)

            SS2=SS2-PSD*(ET-TS)**2
            SS1=SS1-(YL+1.0D0)*DPS(I,J)+PSD*ET
           ENDIF
         ENDDO
      ENDDO

      DQ1=CS*SS1
      DQ2=CS*SS2

      DO J=1,NY
         DO I=1,NX
            IF(PHS(I,J)+XL*DPHS(I,J).NE.0.0D0) THEN
            BB=DPHS(I,J)/(PHS(I,J)+XL*DPHS(I,J))
            DQ1=DQ1+PHT(I,J)*BB
            DQ2=DQ2-PHT(I,J)*BB**2
            ENDIF
         ENDDO
      ENDDO

      XL=XL-DQ1/DQ2
      XL=DMIN1(XLU,XL)

      QQ=HH+CS*SS

C Write out some diagnostics
      IF(VERBOSE) THEN
       WRITE(CHARS,'(''---Iter: '',I3,'//
     :  ''' Accel: '',G15.5,'' Objective function: '',D20.12)')
     :         NIT,XL,QQ
       CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF

C Check for acceleration going negative - if so reset to the 'safe'
C value
      IF(XL.LT.0.0D0) THEN
         CALL UMSPUT('! Warning, acceleration calculated is negative',
     :               1,0,ISTAT)
         XL=XLU
         IF(XL.GT.10.0D0) XL=1.0D0
      ELSE
          IF(DABS(XL-XLP).GT.0.01D0) GO TO 45
      ENDIF

      RETURN
      END

      SUBROUTINE MAXF(PS,DPS,NX,NY,XMAX)
C
C Calculate the maximum acceleration factor obtainable within
C the non-negativity constraint.
C
      INTEGER NX,NY
      DOUBLE PRECISION PS(NX,NY),DPS(NX,NY),XMAX,RT

      INTEGER I,J

      XMAX=1.0D30

      DO J=1,NY
         DO I=1,NX
            IF(DPS(I,J).EQ.0.0D0) GO TO 30
            RT=-PS(I,J)/DPS(I,J)
            IF(RT.LE.0.0D0) GO TO 30
            IF(RT.LT.XMAX) XMAX=RT

30          CONTINUE
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE CIDPSS(PS,PSM,CF,CH,NX,NY,CS,
     :                  SMB,SS,NSTARS,DPSS)
C
C Fill up the DPSS array.
C
C This is for the stationary prior case and is only needed
C if calculating the acceleration factors.
C
      IMPLICIT NONE

      INTEGER NX,NY,I,NSTARS,J
      DOUBLE PRECISION PS(NX,NY),DPSS(NX,NY),SMB,SS,VH,VS
      DOUBLE PRECISION CF(NX,NY),CH(NX,NY)
      DOUBLE PRECISION PSM(NX,NY),CS

C We handle the case of no stars separately
      IF(NSTARS.GT.0) THEN
       DO J=1,NY
          DO I=1,NX
             VH=PS(I,J)*(CF(I,J)-1.0D0)
             VS=CS*(-PS(I,J)*SS-CH(I,J))
     :         +PSM(I,J)*(CF(I,J)-1.0D0)
             DPSS(I,J)=VH+VS
          ENDDO
       ENDDO
      ELSE
       DO J=1,NY
          DO I=1,NX
             VH=PS(I,J)*(CF(I,J)-1.0D0)
             VS=CS*(-PS(I,J)*SS-CH(I,J))
             DPSS(I,J)=VH+VS
          ENDDO
       ENDDO
      ENDIF

      END

      SUBROUTINE FIDPSS(PS,PSM,CF,CFP,CH,NX,NY,CS,
     :                  SS,NSTARS,DPSS)
C
C Fill up the DPSS array.
C
C This is for the floating prior case and is only needed
C if calculating the acceleration factors.
C
C Belatedly modified for the new normalisation, Dec 94
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J,NSTARS
      DOUBLE PRECISION PS(NX,NY),DPSS(NX,NY),SS,VH,VS,VM
      DOUBLE PRECISION CF(NX,NY),CH(NX,NY),CFP(NX,NY)
      DOUBLE PRECISION PSM(NX,NY),CS

C We handle the case of no stars separately
      IF(NSTARS.GT.0) THEN
       DO J=1,NY
          DO I=1,NX
             VH=PS(I,J)*(CF(I,J)-1.0D0)
             VS=CS*PS(I,J)*(CFP(I,J)-1.0D0)+
     :          CS*(-PS(I,J)*SS-CH(I,J))
             VM=PSM(I,J)*(CF(I,J)-1.0D0-CS*SS)
             DPSS(I,J)=VH+VS+VM
          ENDDO
       ENDDO
      ELSE
       DO J=1,NY
          DO I=1,NX
             VH=PS(I,J)*(CF(I,J)-1.0D0)
             VS=CS*PS(I,J)*(CFP(I,J)-1.0D0)+
     :          CS*(-PS(I,J)*SS-CH(I,J))
             DPSS(I,J)=VH+VS
          ENDDO
       ENDDO
      ENDIF

      END
  
      SUBROUTINE FILDPS(DPSS,PSM,CF,NX,NY,CS,SS,DPS)
C
C Fill up the DPS array from the DPSS one.
C
C This is only needed
C if calculating the acceleration factors.
C
C Also changed for the new normalisation, Jan 95
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      DOUBLE PRECISION DPS(NX,NY),PSM(NX,NY)
      DOUBLE PRECISION CF(NX,NY),DPSS(NX,NY),CS,SS

      DO J=1,NY
         DO I=1,NX
            DPS(I,J)=DPSS(I,J)-PSM(I,J)*(CF(I,J)-1.0D0-CS*SS)
         ENDDO
      ENDDO

      END


      SUBROUTINE ACCL(PHT,PS,DPS,PHS,DPHS,NX,NY,CS,XL,
     :                VERBOSE)
C
C Calculate an optimum acceleration factor for point+smooth
C background restorations.
C
C Based on a routine by Leon Lucy.
C
      IMPLICIT NONE

      INTEGER NX,NY

      DOUBLE PRECISION PS(NX,NY),DPS(NX,NY)
      DOUBLE PRECISION PHS(NX,NY),DPHS(NX,NY)
      DOUBLE PRECISION PHT(NX,NY)
      DOUBLE PRECISION CS,XL,XLM1,XLM2,SS,HH,SGP,DSG
      DOUBLE PRECISION XLP,XLU,DQ1,DQ2,QQ,XLM,PSD,PHSD
      DOUBLE PRECISION SS1,SS2,BB

      INTEGER ISTAT,I,J,NIT
      CHARACTER*80 CHARS
      LOGICAL VERBOSE

C First calculate the maximum possible acceleration within the
C positivity restraint
      CALL MAXF(PS,DPS,NX,NY,XLM1)

C and now the same thing in the image plane
      CALL MAXF(PHS,DPHS,NX,NY,XLM2)

      XL=1.0D0

C We now need a 'safe' value for the acceleration - less than
C the max permitted but still positive...
      XLM=DMIN1(XLM1,XLM2)

      IF(VERBOSE) THEN
       WRITE(CHARS,
     :     '(''--Maximum acceleration within non-neg: '','//
     :     'F12.4)') XLM
       CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF

      IF(XLM.GT.1.0D0) THEN
         XLU=1.0D0+0.9D0*(XLM-1.0D0)
      ELSE
         XLU=0.9*XLM
      ENDIF

C Here we do another check - if this value is really very big
C it is likely that we are very close to convergence and using it
C directly will be more harm than good. So we rather arbitrarily
C default to 1.0
      IF(XLU.GT.100.0D0) XLU=1.0D0

      NIT=0

C Return here when iterating
45    CONTINUE

      NIT=NIT+1
      XLP=XL

C Derive the first and second derivatives w.r.t. the
C acceleration factor
      SS=0.0D0
      HH=0.0D0
      SGP=0.0D0
      DSG=0.0D0

      DO J=1,NY
         DO I=1,NX
            SGP=SGP+PS(I,J)+XL*DPS(I,J)
            DSG=DSG+DPS(I,J)
         ENDDO
      ENDDO

      DO J=1,NY
         DO I=1,NX
            PSD=PS(I,J)+XL*DPS(I,J)
            PHSD=PHS(I,J)+XL*DPHS(I,J)

            IF(PSD.GT.0.0D0 .AND. SGP.NE.0.0D0) 
     :            SS=SS-PSD/SGP*DLOG(PSD/SGP)
            IF(PHSD.GT.0.0D0) HH=HH+PHT(I,J)*DLOG(PHSD)
         ENDDO
      ENDDO

      SS1=-SS*DSG/SGP

      DO J=1,NY
         DO I=1,NX
            PSD=PS(I,J)+XL*DPS(I,J)
            IF(PSD.GT.0.0D0 .AND. SGP.NE.0.0D0) 
     :         SS1=SS1-DPS(I,J)*DLOG(PSD/SGP)/SGP
         ENDDO
      ENDDO
             
      SS2=-2.0D0*DSG/SGP*SS1+DSG**2/SGP**2

      DO J=1,NY
         DO I=1,NX
            PSD=PS(I,J)+XL*DPS(I,J)
            IF(PSD.NE.0.0D0) SS2=SS2-(DPS(I,J)/SGP)**2/
     :         ((PS(I,J)+XL*DPS(I,J))/SGP)
         ENDDO
      ENDDO

      DQ1=CS*SS1
      DQ2=CS*SS2

      DO J=1,NY
         DO I=1,NX
            BB=DPHS(I,J)/(PHS(I,J)+XL*DPHS(I,J))
            DQ1=DQ1+PHT(I,J)*BB
            DQ2=DQ2-PHT(I,J)*BB**2
         ENDDO
      ENDDO

      XL=XL-DQ1/DQ2
      XL=DMIN1(XLU,XL)

      QQ=HH+CS*SS

C Write out some diagnostics
      IF(VERBOSE) THEN
       WRITE(CHARS,'(''---Iter: '',I3,'//
     :     ''' Accel: '',G15.5,'' Objective function: '','//
     :     'D20.12)')
     :         NIT,XL,QQ
       CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF

C Check for acceleration going negative - if so reset to the 'safe'
C value
      IF(XL.LT.0.0D0) THEN
         CALL UMSPUT('! Warning, acceleration calculated is negative',
     :               1,0,ISTAT)
         XL=XLU
      ELSE 
         IF(DABS(XL-XLP).GT.0.01D0) GO TO 45
      ENDIF

      RETURN
      END

      SUBROUTINE COPYIM(IN,NX,NY,OUT)
C
C Copy an image to another of the same size
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION IN(NX,NY),OUT(NX,NY)

      INTEGER I,J

      DO J=1,NY
         DO I=1,NX
            OUT(I,J)=IN(I,J)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE INCPSF(RHO,CG,TF,DRHO,CT,TT,NX,NY)
C
C Increment for the PSF when it is being regularized
C From Leon Lucy Code, August 1994
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION RHO(NX,NY),CG(NX,NY),TF(NX,NY)
      DOUBLE PRECISION DRHO(NX,NY)
      DOUBLE PRECISION CT,TT

      DOUBLE PRECISION DRHH,DRHT,DRT,XLM,XL
      CHARACTER*80 CHARS
      INTEGER I,J,ISTAT

      DO J=1,NY
         DO I=1,NX
            DRHH=RHO(I,J)*(CG(I,J)-1.0D0)
            DRT=-TF(I,J)-1.0D0
            DRHT=CT*RHO(I,J)*(DRT-TT+1.0D0)
            DRHO(I,J)=DRHH+DRHT
         ENDDO
      ENDDO

C Find the maximum possible acceleration in the usual way
      CALL MAXF(RHO,DRHO,NX,NY,XLM)

C Display the maximum value
      WRITE(CHARS,'(''-Max. PSF acc.: '',G10.3)') XLM
      CALL UMSPUT(CHARS,1,0,ISTAT)

C Work out a sensible and safe value
      XL=DMIN1(0.9D0*XLM,1.0D0)

      WRITE(CHARS,'(''-Using PSF acc.: '',G10.3)') XL
      CALL UMSPUT(CHARS,1,0,ISTAT)

C and apply it
      DO J=1,NY
         DO I=1,NX
            RHO(I,J)=RHO(I,J)+XL*DRHO(I,J)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SMOOTH(IN,OUT,NX,NY,SIG)
C
C Simple smoothing using direct convolution and a box size
C of 7 by 7 - ie only usable for small SIG
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      DOUBLE PRECISION IN(NX,NY),OUT(NX,NY),SIG,BOX(7,7)
      DOUBLE PRECISION X,Y,T

C First fill the box
      T=0.0D0
      DO J=1,7
         Y=DBLE(J-3)
         DO I=1,7
            X=DBLE(I-3)
            BOX(I,J)=EXP(-(X*X+Y*Y)/(2.0D0*SIG*SIG))
            T=T+BOX(I,J)
         ENDDO
      ENDDO

C Normalise
      DO J=1,7
         DO I=1,7
            BOX(I,J)=BOX(I,J)/T
         ENDDO
      ENDDO

C Convolve
      CALL DIRCON(IN,NX,NY,BOX,7,7,OUT,1)

      RETURN
      END

      SUBROUTINE CHKSUB(DATA,NX,NY,NBX,NBY,ISTAT)
C
C Check the the data sub-sampling is correct - we just look
C at the first few pixels
C
      IMPLICIT NONE

      INTEGER NX,NY,NBX,NBY,I,J,ISTAT
      DOUBLE PRECISION DATA(NX,NY),V

      ISTAT=0

      V=DATA(1,1)
      DO J=1,NBY
         DO I=1,NBX
            IF(DATA(I,J).NE.V) THEN
               ISTAT=1
               RETURN
            ENDIF
         ENDDO
      ENDDO
    
      RETURN 
      END
