C 	Copyright restrictions apply - see st4gem$copyright.stsdas 
C 
      SUBROUTINE ACOADD   !I
CM    PROGRAM MCOADD      !M
C++
C
C ACOADD.F Version 1.6
C
C Lucy co-addition, FFT version with acceleration option.
C
C This code is a preliminary and UNSUPPORTED implementation of
C an algorithm for combining images which have different PSFs.
C
C More details of the method are given in the ST-ECF Newsletter,
C 17, Feb 92, p10 and references given therein.
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
C x_fcoadd.x - the one line of SPP needed for linking purposes.
C
C timtem.x - SPP code for image extension defaulting.
C
C Test version, Richard Hook, ST-ECF, Dec 1991
C               e-mail: rhook@eso.org (Internet)
C                       eso::rhook    (SPAN)
C                       rhook@dgaeso51.bitnet
C
C Note:
C
C This version may also be compiled and linked to run under MIDAS
C as well as IRAF. Essentially all that is required is for the
C lines ending with !I are to be commented out and the comment
C characters removed from lines starting CM. The result should
C then be linked with the additional library 'midint.f' which
C emulates (partially) the F77/VOS routines using MIDAS.
C
C History:
C
C Acceleration added, Richard Hook, ST-ECF, Mar 1992
C Double precision (partly) version, Mar 1992
C Starting image option added, Mar 1992
C Fully double precision version, May 1992
C Weighting error corrected, May 1992
C Additional checks added, May 1992
C Added residual information, June 1992
C Version 1.0, 2nd June 1992
C One D support, 4th June 1992
C Sub-sampling support added, 29th July 1992
C MIDAS compatible version, 7th August 1992
C Version 1.1 frozen, 16th September 1992
C Modified to use faster FFTs and remove powers of 2 restriction, 16/10/92
C Code made more portable (uppercase etc), 21/10/92
C Changed normalisation of coadded output images, 1/12/92
C Increased MAXIM to 99, 4/12/92
C Missing data test additions, 2/2/93 (for version 1.3)
C Simplify parameter file and defaulting, 5/3/93
C Change parameter names to lowercase and tidy up for V1.4, 6/10/93
C Bug fix for DIMS dimensions, V1.5, 14/7/95
C Small change to long character strings which overflow lines to
C avoid warning compilation messages, V1.6, 4/3/98
C--
      IMPLICIT NONE
                   
      INTEGER MAXIMS,USEOF
      PARAMETER (MAXIMS=99)
      PARAMETER (USEOF=-2)
                          
C Iraf global storage
      DOUBLE PRECISION MEMD(1) !I
      COMMON /MEM/MEMD         !I
                                 
C MIDAS global storage
CM    INTEGER MEMD(1)          !M
CM    COMMON /VMR/MEMD         !M
                                 
C Local variables
      INTEGER DATTYP(MAXIMS),DN(MAXIMS),MTYP,NDIMS
      INTEGER DATPNT(MAXIMS),PSFPNT(MAXIMS),MPNT
      INTEGER ISTAT,NX,NY,N,M,NITER,IMSTAT,DIMS(7)
      INTEGER PSFFFT(MAXIMS),WORK,REST,W,DPS,WN,NNEG
      INTEGER CORRFAC,PH(MAXIMS)
      INTEGER IDD(MAXIMS),IDP(MAXIMS),IDDC,IDCOAD(MAXIMS),FID
      DOUBLE PRECISION WEIGHT(MAXIMS),T,FMAX,PT,SCRATCH(4096)
      DOUBLE PRECISION RMSRES,RESMAX
      INTEGER IRMAX,JRMAX
      CHARACTER*80 DATA(MAXIMS),PSF(MAXIMS),MIMAGE,FIM
      CHARACTER*80 COADDR,COADD(MAXIMS),DECON
      CHARACTER*132 CHARS
      CHARACTER*80 IMAGES,PSFS
      INTEGER NIM,NPSF
      INTEGER IMLD,PSFLD,IST,IEND
      INTEGER INR,I,L
      INTEGER NBX,NBY
      DOUBLE PRECISION DRV1,DRV2,DD1,DD2,HH,HHT
      DOUBLE PRECISION XLP,XL1,XL
      LOGICAL ACC
      LOGICAL VERBOSE
      LOGICAL CO
      LOGICAL FIRST
      LOGICAL SUBSAM
      LOGICAL MASK
      LOGICAL RESOUT
C++
   
C Start of code
               
C Midas initialisation
CM    CALL STSPRO('MCOADD') !M
                              
C First get the boolean control logicals
      CALL UCLGSB('verbose',VERBOSE,ISTAT)
      CALL UCLGSB('accel',ACC,ISTAT)
                                    
C Announce the version
      IF(VERBOSE) CALL UMSPUT('+ ACOADD Version 1.6 (March 98)', !I
     :            1,0,ISTAT)                                      !I
CM    IF(VERBOSE) CALL UMSPUT('+ MCOADD Version 1.6 (March 98)', !M
CM   :            1,0,ISTAT)                                      !M
                            
C Initialise the template filename processing
      CALL UCLGST('images',IMAGES,ISTAT)
      CALL TIMOTP(IMAGES,IMLD,ISTAT)
                                    
C Get the names of the input data images
      NIM=1
      IMSTAT=0
      DO I=1,MAXIMS+1
                               
C Check that there are not too many images
         IF(NIM.GT.MAXIMS) THEN
            CALL UMSPUT('! Too many images specified',
     :                  1,0,ISTAT)
            GO TO 99
         ENDIF
              
         CALL TIMXTP(IMLD,DATA(NIM),IMSTAT)
         IF(IMSTAT.EQ.USEOF) GO TO 77
         CALL UIMOPN(DATA(NIM),1,IDD(NIM),ISTAT)
         IF(ISTAT.NE.0) THEN
               CALL UMSPUT('! Unable to open data image',
     :                    1,0,ISTAT)
            GO TO 99
         ELSE
             IF(VERBOSE) CALL UMSPUT
     :         ('-Input Image: '//DATA(NIM),1,0,ISTAT)
         ENDIF
         NIM=NIM+1
      ENDDO
           
 77   CONTINUE
      NIM=NIM-1
               
      CALL TIMCTP(IMLD,ISTAT)
                             
C Open the list of PSF names
      CALL UCLGST('psfs',PSFS,ISTAT)
      CALL TIMOTP(PSFS,PSFLD,ISTAT)
                                   
      NPSF=1
      IMSTAT=0
      DO I=1,MAXIMS+1
         CALL TIMXTP(PSFLD,PSF(NPSF),IMSTAT)
         IF(IMSTAT.EQ.USEOF) GO TO 78
         CALL UIMOPN(PSF(NPSF),1,IDP(NPSF),ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT('! Unable to open PSF image',
     :                    1,0,ISTAT)
            GO TO 99
         ELSE
             IF(VERBOSE) CALL UMSPUT('-Input PSF: '//PSF(NPSF),
     :                   1,0,ISTAT)
         ENDIF
         NPSF=NPSF+1
      ENDDO
           
 78   CONTINUE
      NPSF=NPSF-1
                 
      CALL TIMCTP(PSFLD,ISTAT)
                              
C Do some checks on PSFs size, number etc
      IF(NPSF.NE.NIM) THEN
         CALL UMSPUT('! Different number of PSFs and Images',
     :                 1,0,ISTAT)
         GO TO 99
      ELSE
         WRITE(CHARS,'(''-There are '',I4,'' images/PSFs'')') NIM
         IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
                                                 
C In either case we set the acceleration to 1.0 as default
         XL=1.0D0
      ENDIF
           
C Get the shapes and sizes of the images
      CALL UIMGID(IDD(1),DATTYP(1),NDIMS,DIMS,ISTAT)
                                                    
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

      DO N=1,NIM
        CALL UIMGID(IDD(N),DATTYP(N),NDIMS,DIMS,ISTAT)
                                                      
         IF(NDIMS.EQ.1) THEN
            DIMS(2)=1
         ENDIF
              
C Check they are OK
        IF(DIMS(1).NE.NX .OR.
     :     DIMS(2).NE.NY) THEN
           CALL UMSPUT('! Data/PSF size/shape error',1,0,ISTAT)
           WRITE(CHARS,'('' N,NX,NY,DIMS(1),DIMS(2): '',5I6)')
     :        N,NX,NY,DIMS(1),DIMS(2)
           CALL UMSPUT(CHARS,1,0,ISTAT)
           GO TO 99
        ENDIF
      ENDDO
           
C Same for the PSFs
      DO N=1,NIM
        CALL UIMGID(IDP(N),DATTYP(N),NDIMS,DIMS,ISTAT)
                                                      
         IF(NDIMS.EQ.1) THEN
            DIMS(2)=1
         ENDIF
              
C Check they are OK
        IF(NX.NE.DIMS(1) .OR.
     :     NY.NE.DIMS(2)) THEN
           CALL UMSPUT('! Data/PSF size/shape error',1,0,ISTAT)
           WRITE(CHARS,'('' N,NX,NY,DIMS(1),DIMS(2): '',5I6)')
     :        N,NX,NY,DIMS(1),DIMS(2)
           CALL UMSPUT(CHARS,1,0,ISTAT)
           GO TO 99
        ENDIF
      ENDDO
           
C Check for the existence of a mask image, if so open it,
C if not just set the flag
      CALL UCLGST('maskim',MIMAGE,ISTAT)
      IF(MIMAGE.EQ.' ') THEN
         MASK=.FALSE.
      ELSE
        
         CALL UIMOPN(MIMAGE,1,MPNT,ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT('! Unable to open mask image',1,0,ISTAT)
           GO TO 99
         ENDIF
         MASK=.TRUE.

         IF(VERBOSE) CALL UMSPUT('-Mask Image: '//MIMAGE,
     :                           1,0,ISTAT)

C Check that the mask is the right size/shape
         CALL UIMGID(MPNT,MTYP,NDIMS,DIMS,ISTAT)
         IF(NDIMS.EQ.1) DIMS(2)=1
      
         IF(NX.NE.DIMS(1) .OR.
     :      NY.NE.DIMS(2)) THEN
            CALL UMSPUT('! Mask size/shape error',1,0,ISTAT)
            GO TO 99
         ENDIF

C Check that we aren't also trying to use acceleration - this is not
C currently implemented
         IF(ACC) THEN
            CALL UMSPUT('! The accelerated algorithm may not be'//
     :                  ' used with a mask',1,0,istat)
            GO TO 99
         ENDIF
      ENDIF

C Allocate space for the data
      DO N=1,NIM
        CALL UDMGET(NX*NY,7,DATPNT(N),ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for data array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
             
C Read the data into the memory
        CALL UIGS2D(IDD(N),1,NX,1,NY,MEMD(DATPNT(N)),ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to read in data input image',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
             
C Check that the data is non-negative
        CALL ZAPNEG(MEMD(DATPNT(N)),NX,NY,NNEG)
        IF(NNEG.GT.0) THEN
             WRITE(CHARS,
     :       '(''! Warning, data image #'',I2,'//
     :       ''' contains'',I7,'' neg. values, these'','//
     :       ''' have been zeroed.'')') N,NNEG
            CALL UMSPUT(CHARS,1,0,ISTAT)
        ENDIF
             
C Allocate space for the PSFs
        CALL UDMGET(NX*NY,7,PSFPNT(N),ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for PSF array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
             
C Read the PSF into the memory
        CALL UIGS2D(IDP(N),1,NX,1,NY,MEMD(PSFPNT(N)),ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to read in PSF input image',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF

C Get working space arrays for the FFTs
        CALL UDMGET(NX*NY*2,7,PSFFFT(N),ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for PSF FFT array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
      ENDDO
             
C Allocate space for the mask, if there is one
      IF(MASK) THEN
        CALL UDMGET(NX*NY,7,WN,ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for mask array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF

C Read the mask into the memory
        CALL UIGS2D(MPNT,1,NX,1,NY,MEMD(WN),ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to read in mask image',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
      ENDIF

C Create the output data arrays
      CALL UCLGST('decon',DECON,ISTAT)
      IF(DECON.EQ.' ') THEN
         RESOUT=.FALSE.
      ELSE
         RESOUT=.TRUE.
         CALL UIMCRE(DECON,7,2,DIMS,IDDC,ISTAT)
         IF(ISTAT.NE.0) THEN
            CALL UMSPUT(
     :        '! Unable to open output deconvolved image',
     :                    1,0,ISTAT)
            GO TO 99
         ELSE
            IF(VERBOSE) CALL UMSPUT(
     :        '-Created output (deconvolved) image: '//DECON,
     :        1,0,ISTAT)
         ENDIF
      ENDIF
           
C Get the root name for the coadded output images
      CALL UCLGST('coaddr',COADDR,ISTAT)

C Check if a name is given
      IF(COADDR.EQ.' ') THEN
         CO=.FALSE.

C Check if there are no output files at all - if so abort
         IF(.NOT.RESOUT) THEN
          CALL UMSPUT('! No output image names given - aborting',
     :                  1,0,ISTAT)
          GO TO 99
         ENDIF
      ELSE

         CO=.TRUE.
                                        
C Find the length of the string to avoid blanks in the names
         CALL LENSTR(COADDR,IST,IEND)
         L=IEND-IST+1
                                  
         DO N=1,NIM
            COADD(N)(1:L)=COADDR(IST:IEND)
            COADD(N)(L+1:L+1)='_'
            IF(N.LT.10) THEN
               WRITE(COADD(N)(L+2:L+2),'(I1)') N
            ELSE
               WRITE(COADD(N)(L+2:L+3),'(I2)') N
            ENDIF
            CALL UIMCRE(COADD(N),7,2,DIMS,IDCOAD(N),ISTAT)
            IF(ISTAT.NE.0) THEN
               CALL UMSPUT('! Unable to open output coadded image',
     :                       1,0,ISTAT)
               GO TO 99
            ELSE
             IF(VERBOSE) CALL UMSPUT(
     :           '-Created output (coadded) image: '//COADD(N),
     :           1,0,ISTAT)
            ENDIF
          ENDDO
      ENDIF
            
89    CONTINUE
              
C Get the sub-sampling values in X and Y and check that
C the X and Y dimensions are multiples of them.
      CALL UCLGSI('xsubsam',NBX,ISTAT)
      CALL UCLGSI('ysubsam',NBY,ISTAT)
                                      
      IF(MOD(NX,NBX).NE.0 .OR.
     :   MOD(NY,NBY).NE.0) THEN
         CALL UMSPUT('! Warning - arrays are not multiples'//
     :               ' of the sub-sampling',1,0,ISTAT)
      ENDIF
           
      IF(NBX.NE.1 .OR. NBY.NE.1) THEN
         SUBSAM=.TRUE.

C Check that we aren't also using a mask - this is currently
C not implemented
         IF(MASK) THEN
            CALL UMSPUT('! Subsampling may not be used with'//
     :                  ' a mask',1,0,ISTAT)
            GO TO 99
         ENDIF
      ELSE
         SUBSAM=.FALSE.
      ENDIF
           
C Prepare the arrays for the FFT routine
C and do the FFTs of the PSF and rotated PSF
C
C At this point we also check that the PSFs are normalised to
C a total of 1 - if they aren't we warn the user but continue
      DO N=1,NIM
                
C First check for negative values, if there are any set them
C to zero
         CALL ZAPNEG(MEMD(PSFPNT(N)),NX,NY,NNEG)
         IF(NNEG.GT.0) THEN
             WRITE(CHARS,
     :       '(''! Warning, PSF #'',I2,'//
     :       ''' contained '',I7,'' neg. values, these'','//
     :       ''' have been set to zero.'')') N,NNEG
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF
              
         CALL TOTAL(MEMD(PSFPNT(N)),NX,NY,PT)
         IF(DABS(PT-1.0D0).GT.1.0E-5) THEN
            WRITE(CHARS,
     :       '(''! Warning, PSF #'',I2,'//
     :       ''' not normalised to a total of 1.0'','//
     :       '''. Sum is '',G15.5)') N,PT
            CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF
              
         CALL DFILL(MEMD(PSFPNT(N)),NX,NY,MEMD(PSFFFT(N)))
         CALL DFOURT(MEMD(PSFFFT(N)),DIMS,2,-1,0,SCRATCH)
                                                         
C At this point we have finished with the PSFs so they can be closed
C and the memory freed.
         CALL UIMCLO(IDP(N),ISTAT)
         CALL UDMFRE(PSFPNT(N),7,ISTAT)
      ENDDO
           
C Create other arrays needed, workspace etc
        CALL UDMGET(NX*NY,7,W,ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for working array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
             
        CALL UDMGET(NX*NY*2,7,WORK,ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for workspace array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
             
        CALL UDMGET(NX*NY,7,REST,ISTAT)
        IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for restored array',
     :                 1,0,ISTAT)
           GO TO 99
        ENDIF
             
C PHI arrays and DPS are only needed in the accelerated case
       IF(ACC)THEN
        DO N=1,NIM
           CALL UDMGET(NX*NY,7,PH(N),ISTAT)
           IF(ISTAT.NE.0) THEN
              CALL UMSPUT(
     :           '! Unable to allocate memory for Phi array',
     :                    1,0,ISTAT)
              GO TO 99
           ENDIF
         ENDDO
              
         CALL UDMGET(NX*NY,7,DPS,ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to allocate memory for delta Psi array',
     :                 1,0,ISTAT)
           GO TO 99
         ENDIF
       ENDIF
            
C Get space for and initialise the correction factor array
       CALL UDMGET(NX*NY,7,CORRFAC,ISTAT)
       IF(ISTAT.NE.0) THEN
          CALL UMSPUT('! Unable to create corr. fac. image',
     :                    1,0,ISTAT)
          GO TO 99
       ENDIF
            
C Allocate space for the DN arrays if there is a mask
      IF(MASK) THEN
         DO M=1,NIM
            CALL UDMGET(NX*NY,7,DN(M),ISTAT)
            IF(ISTAT.NE.0) THEN
               CALL UMSPUT(
     :         '! Unable to allocate memory for DN array',
     :                    1,0,ISTAT)
               GO TO 99
            ENDIF

C We also need to multiply the PHT arrays by the mask
            CALL MULT(MEMD(WN),MEMD(DATPNT(M)),NX,NY,
     :                MEMD(DATPNT(M)))
         ENDDO
      ENDIF

C Calculate the weights - just the normalised total counts in each image
       T=0.0D0
       DO N=1,NIM
          CALL TOTAL(MEMD(DATPNT(N)),NX,NY,WEIGHT(N))
          T=T+WEIGHT(N)
       ENDDO
            
       DO N=1,NIM
          WRITE(CHARS,'(''-Image '',I2,'' Total counts: '','//
     :       'G10.4,'' Weight: '',F5.3)') N,WEIGHT(N),WEIGHT(N)/T
          IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
          WEIGHT(N)=WEIGHT(N)/T
       ENDDO
            
C Check whether we have a first estimate image, if so try to get it
      CALL UCLGST('firstim',FIM,ISTAT)
                                         
      IF(FIM.EQ.' ') THEN
         FIRST=.FALSE.
          
C Prepare the first estimate - just a flat image with the same total flux
          CALL FILCON(MEMD(REST),NX,NY,T/DFLOAT(NX*NY))
      ELSE
         CALL UIMOPN(FIM,1,FID,ISTAT)
         IF(ISTAT.NE.0) THEN
           CALL UMSPUT(
     :        '! Unable to open first estimate image file',
     :                 1,0,ISTAT)
           GO TO 99
         ELSE
            IF(VERBOSE) CALL UMSPUT
     :         ('-First Estimate Image: '//FIM,1,0,ISTAT)
         ENDIF
              
C Read it in and then close the image
         CALL UIGS2D(FID,1,NX,1,NY,MEMD(REST),ISTAT)
         CALL UIMCLO(FID,ISTAT)
      ENDIF
           
C If there is a mask we can prepare the DN arrays by correlating the
C PSF with the mask
      IF(MASK) THEN
         DO M=1,NIM
            CALL DCONV(MEMD(WN),NX,NY,MEMD(WORK),
     :                    MEMD(PSFFFT(M)),MEMD(DN(M)),-1)
         ENDDO
      ENDIF

C Now we can start the iterative procedure
      CALL UCLGSI('niter',NITER,ISTAT)
                                      
      DO N=1,NITER
                  
         WRITE(CHARS,'(''--Starting iteration '',I5)') N
         IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
                                                 
C The correction array has to be initialised each time
        CALL FILCON(MEMD(CORRFAC),NX,NY,0.0D0)
                                              
        DO M=1,NIM
                  
C We must first multiply the current estimate by the weight
C to get something with the same total flux as the data 
         CALL MULC(MEMD(REST),NX,NY,WEIGHT(M),MEMD(W))
                                                      
C Convolve the current estimate with the PSF to give Phi
         IF(ACC) THEN
          CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :              MEMD(PSFFFT(M)),MEMD(PH(M)),1)
         ELSE
          CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :              MEMD(PSFFFT(M)),MEMD(W),1)
         ENDIF
              
C Calculate and display the residual information to assist closeness
C of fit assessment
         IF(VERBOSE) THEN
          IF(MASK) THEN
           IF(ACC) THEN
             CALL RESINM(MEMD(PH(M)),MEMD(DATPNT(M)),
     :             MEMD(WN),NX,NY,RMSRES,RESMAX,IRMAX,JRMAX)
           ELSE
             CALL RESINM(MEMD(W),MEMD(DATPNT(M)),
     :             MEMD(WN),NX,NY,RMSRES,RESMAX,IRMAX,JRMAX)
           ENDIF
          ELSE
           IF(ACC) THEN
             CALL RESINF(MEMD(PH(M)),MEMD(DATPNT(M)),
     :                    NX,NY,RMSRES,RESMAX,IRMAX,JRMAX)
           ELSE
             CALL RESINF(MEMD(W),MEMD(DATPNT(M)),
     :                    NX,NY,RMSRES,RESMAX,IRMAX,JRMAX)
           ENDIF
          ENDIF
               
          WRITE(CHARS,
     : '(''--Image #'',I2,'' RMS residual: '',D10.4,'//
     :                 ''' Max residual of '',D10.4,'//
     :                 ''' at ('',I4'','',I4,'')'')')
     :             M,RMSRES,RESMAX,IRMAX,JRMAX
          CALL UMSPUT(CHARS,1,0,ISTAT)
         ENDIF
              
C If we are sub-sampling we must do so here.
      IF(SUBSAM) THEN
         IF(ACC) THEN
            CALL REBIN(MEMD(PH(M)),NX,NY,NBX,NBY)
         ELSE
            CALL REBIN(MEMD(W),NX,NY,NBX,NBY)
         ENDIF
      ENDIF
C Divide the data by this convolution
         IF(ACC) THEN
          CALL DIVIDE(MEMD(DATPNT(M)),MEMD(PH(M)),NX,NY,MEMD(W))
         ELSE
          CALL DIVIDE(MEMD(DATPNT(M)),MEMD(W),NX,NY,MEMD(W))
         ENDIF
              
C Convolve again, this time by the rotated PSF to make the
C correction factor array
         CALL DCONV(MEMD(W),NX,NY,MEMD(WORK),
     :             MEMD(PSFFFT(M)),MEMD(W),-1)
                                              
C If there is a mask this must be divided by DN
         IF(MASK) THEN
            CALL DIVIDE(MEMD(W),MEMD(DN(M)),NX,NY,MEMD(W))
         ENDIF

C Add the correction factor so deduced into the overall correction
C factor with the correct weighting
         CALL ADDW(MEMD(W),WEIGHT(M),NX,NY,MEMD(CORRFAC))
        ENDDO
             
C In the accelerated case calculate the maximum possible multiplication
C factor which is consistent with the non-negativity contraint
        IF(ACC) THEN
           CALL FIMAXF(MEMD(REST),MEMD(CORRFAC),NX,NY,FMAX)
           WRITE(CHARS,
     :   '(''--Max. speed-up possible'','//
     :     ''' within non-negativity (fmax): '','//
     :      'F10.5)') FMAX
           IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
                                                   
C Work out the increment in Psi
           CALL PHINC(MEMD(REST),MEMD(CORRFAC),NX,NY,MEMD(DPS))
                                                               
C We now start the iteration loop for the Newton-Raphson search
                                                               
           IF(VERBOSE) THEN
             CALL UMSPUT('---Starting search for optimal speed'//
     :                 '-up factor...',1,0,ISTAT)
           ENDIF
                
           XL1=1.0D0
           INR=0
                
145        CONTINUE
           XLP=XL1
                  
           DRV1=0.0D0
           DRV2=0.0D0
           HHT=0.0D0
                    
           DO M=1,NIM
                     
C Convolve the increment in psi to get the increment in phi
              CALL DCONV(MEMD(DPS),NX,NY,MEMD(WORK),
     :                  MEMD(PSFFFT(M)),MEMD(W),1)
                                                  
C Multiply by the weight again to conserve the flux total
C (this is again a new addition)
              CALL MULC(MEMD(W),NX,NY,WEIGHT(M),MEMD(W))
                                                        
C If we are sub-sampling we must do so here again.
      IF(SUBSAM) THEN
         CALL REBIN(MEMD(W),NX,NY,NBX,NBY)
      ENDIF
           
C Calculate the derivatives needed for Newton-Raphson
              CALL DERIVS(MEMD(DATPNT(M)),MEMD(PH(M)),MEMD(W),
     :                    NX,NY,XL1,DD1,DD2,HH)
                                               
              HHT=HH+HHT
                        
C Add in the derivs with the correct weights
              DRV1=DRV1+DD1*WEIGHT(M)
              DRV2=DRV2+DD2*WEIGHT(M)
           ENDDO
                
           INR=INR+1
                    
C Calculate next estimate
           XL1=XL1-DRV1/DRV2
                            
           WRITE(CHARS,
     :    '(''---Iteration: '',I3,'' Speed up: '',F9.4,'//
     :    ''' Log.Likelihood: '','//
     :    'D15.8)') INR,XL1,HHT
           IF(VERBOSE)  CALL UMSPUT(CHARS,1,0,ISTAT)
                                                    
C Check for negative likelihood gradient
           IF(INR.EQ.1 .AND.DRV1.LT.0.0D0) THEN
              CALL UMSPUT('! Error - drv1<0 for inr=1',1,0,ISTAT)
              CALL UMSPUT('! No acceleration this time round',1,0,
     :                    ISTAT)
              XL1=1.0D0
              GO TO 148
           ENDIF
                
C Check for speed-up too big for non-negativity
C In which case set the speed-up-factor at 70% of the max
C permitted factor
           IF(XL1.GT.FMAX) THEN
              XL1=1.0D0+0.7D0*(FMAX-1.0D0)
              IF(VERBOSE) CALL UMSPUT(
     :        '---Avoiding non-negativity violation'//
     :        ' (accel=1.0+0.7*(fmax-1))',
     :                    1,0,ISTAT)
              GO TO 148
           ENDIF
                
           IF(ABS(XL1-XLP).GT.0.02D0) GO TO 145
                                               
148        CONTINUE
           XL=XL1
           WRITE(CHARS,'(''---Speed up factor used: '','//
     :           'F10.5)') XL
           IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
        ENDIF
             
C Multiply the last estimate by the total correction factor from all the
C images and hence obtain the next estimate
C We also remormalise at this stage to prevent accumulating
C normalisation errors - but only if there isn't a mask - if there
C is we DON'T renormalise
        CALL UPCORR(MEMD(REST),MEMD(CORRFAC),XL,NX,NY,
     :                MASK,MEMD(REST),T,VERBOSE)
      ENDDO
           
C Write out the results, deconvolved and then the coadded images
C which are the convolutions with the relavant PSFs
      IF(RESOUT) THEN
         IF(VERBOSE)
     :      CALL UMSPUT('-Writing output deconvolved image: '//DECON
     :              ,1,0,ISTAT)
                               
          CALL UIPS2D(IDDC,1,NX,1,NY,MEMD(REST),ISTAT)
          CALL UIMCLO(IDDC,ISTAT)
      ENDIF
                              
       IF(CO) THEN
          DO N=1,NIM
             CALL DCONV(MEMD(REST),NX,NY,MEMD(WORK),
     :              MEMD(PSFFFT(N)),MEMD(W),1)
                                              
                                                    
             IF(VERBOSE)
     :        CALL UMSPUT('-Writing output co-added image: '//
     :                    COADD(N),1,0,ISTAT)
                                          
             CALL UIPS2D(IDCOAD(N),1,NX,1,NY,MEMD(W),ISTAT)
             CALL UIMCLO(IDCOAD(N),ISTAT)
          ENDDO
       ENDIF
            
99     CONTINUE
               
C Close all the input data images
      DO N=1,NIM
         CALL UIMCLO(IDD(N),ISTAT)
      ENDDO
           
C Free all the dynamic arrays
      DO N=1,NIM
         CALL UDMFRE(DATPNT(N),7,ISTAT)
         CALL UDMFRE(PSFFFT(N),7,ISTAT)
         IF(ACC) CALL UDMFRE(PH(N),7,ISTAT)
      ENDDO
           
      CALL UDMFRE(CORRFAC,7,ISTAT)
      CALL UDMFRE(WORK,7,ISTAT)
      CALL UDMFRE(REST,7,ISTAT)
      CALL UDMFRE(W,7,ISTAT)
      IF(ACC) CALL UDMFRE(DPS,7,ISTAT)
                                      
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
         
      SUBROUTINE UPCORR(IN1,IN2,XL,NX,NY,MASK,
     :                  OUT,T,VERBOSE)
C
C Apply the update, with a speed up factor
C
C Also renormalise to prevent accumulating errors
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J,ISTAT,NNEG
      DOUBLE PRECISION IN1(NX,NY),IN2(NX,NY)
      DOUBLE PRECISION T,V,OUT(NX,NY),XL,TOTAL
      CHARACTER*80 CHARS
      LOGICAL VERBOSE,MASK
                     
      TOTAL=0.0D0
      NNEG=0
            
C Work out the total counts allowing for negative number
C corrections
      DO J=1,NY
         DO I=1,NX
            V=IN1(I,J)*(1.0D0+XL*(IN2(I,J)-1.0D0))
            TOTAL=TOTAL+V
            IF(V.LT.0.0D0) THEN
               TOTAL=TOTAL-V
               OUT(I,J)=0.0D0
               NNEG=NNEG+1
            ELSE
               OUT(I,J)=V
            ENDIF
         ENDDO
      ENDDO
           
C Do the renormalisation-if there is no mask
      IF(MASK) THEN
         WRITE(CHARS,'(''--Total in restored image: '','//
     :         'G20.5)') TOTAL
         IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
      ELSE
         DO J=1,NY
            DO I=1,NX
               OUT(I,J)=OUT(I,J)*T/TOTAL
            ENDDO
         ENDDO
           
         WRITE(CHARS,
     : '(''--Renormalising restored image (Factor: '',F20.18,'')'')')
     :    T/TOTAL
         IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
                                              
      ENDIF

      IF(NNEG.GT.0) THEN
         WRITE(CHARS,'(''--A total of '',I8,'//
     :         ''' Negative points set to zero.'')') NNEG
         IF(VERBOSE) CALL UMSPUT(CHARS,1,0,ISTAT)
      ENDIF
           
      RETURN
      END

      SUBROUTINE MULT(IN1,IN2,NX,NY,OUT)
C
C Just multiply one array by another of the same size
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      DOUBLE PRECISION IN1(NX,NY),IN2(NX,NY),OUT(NX,NY)

      DO J=1,NY
         DO I=1,NX
            OUT(I,J)=IN1(I,J)*IN2(I,J)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE DIVIDE(IN1,IN2,NX,NY,OUT)
C
C Just divide one array by another of the same size
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J
      DOUBLE PRECISION IN1(NX,NY),IN2(NX,NY),OUT(NX,NY)
                                                       
      DO J=1,NY
         DO I=1,NX
            IF(IN2(I,J).EQ.0.0D0) THEN
               OUT(I,J)=0.0D0
            ELSE
               OUT(I,J)=IN1(I,J)/IN2(I,J)
            ENDIF
         ENDDO
      ENDDO
           
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
         
      SUBROUTINE ADDW(DATA,W,NX,NY,OUT)
C
C Add an image, weighted, to another
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J
      DOUBLE PRECISION DATA(NX,NY),W,OUT(NX,NY)
                                               
      DO J=1,NY
         DO I=1,NX
            OUT(I,J)=DATA(I,J)*W+OUT(I,J)
         ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE FIMAXF(PS,CF,NX,NY,FMAX)
C
C Find the maximum multiplying factor possible for acceleration
C which is consistent with the non-negativity contraint.
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,I,J
      DOUBLE PRECISION PS(NX,NY),CF(NX,NY),FMAX,RT
                                                  
      FMAX=1.0D10
                 
      DO J=1,NY
         DO I=1,NX
             IF(PS(I,J).LT.1.0D-10 .OR.
     :          CF(I,J).EQ.1.0) GO TO 88
             RT=-1.0D0/(CF(I,J)-1.0D0)
             IF(RT.LT.0.0D0) GO TO 88
             IF(RT.LT.FMAX) FMAX=RT
 88          CONTINUE
          ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE DERIVS(PHT,PH,DPH,NX,NY,XL,DD1,DD2,HH)
C
C Calculate derivatives needed in Newton-Raphson search for
C optimal acceleration factor.
C
C Based on code by Leon Lucy
C
      IMPLICIT NONE
                   
      INTEGER NX,NY
      DOUBLE PRECISION PHT(NX,NY),PH(NX,NY),DPH(NX,NY)
      DOUBLE PRECISION XL
      DOUBLE PRECISION DD1,DD2,HH
                                 
      INTEGER I,J
      DOUBLE PRECISION DM,AA,DDPH,DPHT
                                      
      DD1=0.0D0
      DD2=0.0D0
      HH=0.0D0
              
C NB here we are ignoring negative values of dm - these
C may represent pathological situations which should be checked
C out.
      DO J=1,NY
         DO I=1,NX
            DDPH=DPH(I,J)
            DM=PH(I,J)+XL*DDPH
            IF(DM.GT.1.0D-20) THEN
               DPHT=PHT(I,J)
               AA=DDPH/DM
               DD1=DD1+DPHT*AA
               DD2=DD2-DPHT*AA*AA
               HH=HH+DPHT*DLOG(DM)
            ENDIF
         ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE PHINC(PS,CF,NX,NY,DPS)
C
C Calculate the increment in Psi, this is only needed in the
C accelerated case
C
      IMPLICIT NONE
                   
      INTEGER NX,NY
      DOUBLE PRECISION PS(NX,NY),CF(NX,NY),DPS(NX,NY)
      INTEGER I,J
                 
      DO J=1,NY
         DO I=1,NX
            DPS(I,J)=PS(I,J)*(CF(I,J)-1.0D0)
         ENDDO
      ENDDO
           
      RETURN
      END
         
      SUBROUTINE LENSTR(STRING,I1,I2)
C
C Find the start and end of a string
C
      IMPLICIT NONE
                   
      CHARACTER*(*) STRING
      INTEGER I,I1,I2
      LOGICAL IN
                
      IN=.FALSE.
                
      DO I=1,LEN(STRING)
         IF(STRING(I:I).NE.' ' .AND.
     :      .NOT.IN) THEN
           I1=I
           IN=.TRUE.
         ENDIF
              
         IF(STRING(I:I).EQ.' ' .AND.
     :      IN) THEN
            I2=I-1
            GO TO 99
         ENDIF
      ENDDO
           
99    CONTINUE
              
      RETURN
      END
         
      SUBROUTINE RESINM(IM1,IM2,M,NX,NY,RMSRES,RESMAX,
     :                   IRMAX,JRMAX)
C
C Calculate the RMS residual between two images and the
C largest residual. Version for masks.
C
      IMPLICIT NONE

      INTEGER NX,NY,IRMAX,JRMAX
      DOUBLE PRECISION RMSRES,RESMAX,IM1(NX,NY),IM2(NX,NY)
      DOUBLE PRECISION M(NX,NY)

      INTEGER I,J,NV
      DOUBLE PRECISION RES,T

      T=0.0D0
      RESMAX=0.0D0
      NV=0

      DO J=1,NY
         DO I=1,NX
            IF(M(I,J).NE.0) THEN
               NV=NV+1
               RES=IM1(I,J)-IM2(I,J)
               T=T+RES*RES
               IF(DABS(RES).GT.DABS(RESMAX)) THEN
                  RESMAX=RES
                  IRMAX=I
                  JRMAX=J
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      RMSRES=DSQRT(T/DBLE(NV))

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
         
      SUBROUTINE ZAPNEG(DATA,NX,NY,NNEG)
C
C Find negative points in an array and set them to zero
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,NNEG
      DOUBLE PRECISION DATA(NX,NY)
      INTEGER I,J
                 
      NNEG=0
            
      DO J=1,NY
         DO I=1,NX
            IF(DATA(I,J).LT.0.0D0) THEN
               DATA(I,J)=0.0D0
               NNEG=NNEG+1
            ENDIF
         ENDDO
      ENDDO
           
      RETURN
      END
