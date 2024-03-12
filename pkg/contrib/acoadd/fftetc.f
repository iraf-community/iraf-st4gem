C     Copyright restrictions apply - see st4gem$copyright.stsdas 
C 
C FFTETC.F - some useful FFT and convolution routines.
C
C This file contains three subroutines:
C
C  i) An FFT implementation which is fast and flexible
C     (even though it is 25 years old). It is called DFOURT.
C
C ii) A simple piece of code which puts an image into the
C     correct form to be swallowed by DFOURT. It is DFILL.
C
C iii) A convolution code which uses the former to convolve a
C      supplied image with the FFT of a supplied PSF. It is
C      DCONV.
C
C The authors of the first of these are given in the comment
C block at the start of the code. The other two smaller routines
C were written by Richard Hook, Autumn 1992.
C
C Note that all these routines work fully in DOUBLE PRECISION.
C
C Note - this is a stripped down version, the full one
C        with different kinds of FFT and single precisio
C        versions is in 'bigfftetc.f'
C
CCC Note--------
C Note - as of October 1993 a bug was fixed in the DFILL
C routine. This affects the first row and column and normally has
C no effect. "Normally" means cases where the PSF is zero at the
C edge of the frame. However if you have been using an earlier version
C please change to this one.
CCC Note--------
C
      SUBROUTINE DFILL(DATA,NX,NY,OUT)
C
C Fill up an array in the way required by the DFOURT
C FFT routines
C
C Supplied:
C  
C DATA - a double precision 2d array.
C
C NX,NY - integers, the array dimensions.
C
C Returned:
C
C OUT - a double precision 1d array in the correct form
C       to be used by the DFOURT FFT routines. The data values
C       in DATA are simply put into alternating elements and
C       the rest set to zero.
C
C Note this routine has a 'fudge factor' of 1 pixel to align coordinate
C systems correctly.
C
      IMPLICIT NONE
                   
      INTEGER NX,NY,IP,I,J
      DOUBLE PRECISION DATA(NX,NY)
      DOUBLE PRECISION OUT(NX*NY*2)
                                   
      IP=1
          
      IF(NY.GT.1) THEN
         OUT(IP)=DATA(NX,NY)
         OUT(IP+1)=0.0D0
         IP=IP+2
         DO I=1,NX-1
            OUT(IP)=DATA(I,NY)
            OUT(IP+1)=0.0D0   
            IP=IP+2
         ENDDO
              
         DO J=1,NY-1
            OUT(IP)=DATA(NX,J)
            OUT(IP+1)=0.0D0
            IP=IP+2
            DO I=1,NX-1
               OUT(IP)=DATA(I,J)
               OUT(IP+1)=0.0D0
               IP=IP+2
            ENDDO
          ENDDO
       ELSE
         DO I=1,NX
            OUT(IP)=DATA(I,1)
            OUT(IP+1)=0.0D0
            IP=IP+2
         ENDDO
      ENDIF
           
      RETURN
      END
         
      SUBROUTINE DCONV(DATA,NX,NY,WORK,PSFFFT,OUT,FLAG)
C
C Convolve a one or two dimensional REAL array with a PSF which is supplied as
C a double precision FFT. Uses the DFOURT FFT code. 
C
C Supplied:
C
C DATA - 2d double precision image array of dimensions (NX,NY)
C
C NX,NY - integers, the dimensions of this array
C
C WORK - double precision 1d array of dimension (NX*NY*2)
C        This is used as workspace for the FFTs.
C
C PSFFFT - double precision 1d array of dimension (NX*NY*2)
C          This is the FFT of the PSF as created by the routine
C          DFOURT.
C
C FLAG - integer (1 or -1). If this is set to 1 a simple
C        convolution operation will be performed. If -1
C        a correlation will be done instead. This is equivalent
C        to using a PSF rotated by 180 degrees or multiplying the
C        imaginary part of the PSF FFT by -1.
C
C Returned:
C
C OUT - 2d double precision image array of dimensions (NX,NY)
C       The result of the convolution.
C
C Richard Hook, ST-ECF, October 1991
C Double precision version, March 1992
C 1D version, August 1992
C Modified to use DFOURT and improved comments, October 1992
C
      IMPLICIT NONE
                   
      INTEGER I,J,N,IP,IQ
      DOUBLE PRECISION A,B,C,D,FAC
      INTEGER NX,NY                ! ARRAY SIZE
      INTEGER NDIMS,DIMS(2),NX2,NY2
      INTEGER FLAG                 ! ROTATED OR NOT?
      DOUBLE PRECISION DATA(NX,NY)             ! INPUT DATA ARRAY
      DOUBLE PRECISION OUT(NX,NY)
      DOUBLE PRECISION PSFFFT(NX*NY*2)         ! FFT OF PSF
      DOUBLE PRECISION WORK(NX*NY*2)           ! WORKSPACE FOR FFT OF DATA
      DOUBLE PRECISION SCRATCH(40960)
      DOUBLE PRECISION DFLAG
                            
      DFLAG=DBLE(FLAG)
      IF(NY.EQ.1) THEN
         NDIMS=1
      ELSE
         NDIMS=2
      ENDIF
           
C Copy the data into the real part of the working array
      IP=1
      DO J=1,NY
         DO I=1,NX
            WORK(IP)=DATA(I,J)
            WORK(IP+1)=0.0D0
            IP=IP+2
         ENDDO
      ENDDO
           
C Take the forward FFT of the data array
      DIMS(1)=NX
      DIMS(2)=NY
                
      CALL DFOURT(WORK,DIMS,NDIMS,-1,0,SCRATCH)
                                               
C Do a complex multiplication of the FFT of the data and of the PSF
C
C If the rotation flag is set (to -1) we multiply the imaginary part
C of the FFT by -1
      IP=1
      IQ=2
      DO N=1,NX*NY
         A=WORK(IP)
         B=WORK(IQ)
         C=PSFFFT(IP)
         D=DFLAG*PSFFFT(IQ)
         WORK(IP)=A*C-B*D
         WORK(IQ)=A*D+B*C
         IP=IP+2
         IQ=IQ+2
      ENDDO
           
C Take the inverse FFT
      CALL DFOURT(WORK,DIMS,NDIMS,1,1,SCRATCH)
                                              
C Copy the result into the output array
C and 'quadrant swap':
      FAC=DFLOAT(NX*NY)
      IP=1
          
      NX2=MAX(1,NX/2)
      NY2=MAX(1,NY/2)
                     
      IF(NY.GT.1) THEN
       DO J=1,NY2
          DO I=1,NX2
             OUT(I+NX/2,J+NY/2)=WORK(IP)/FAC
             IP=IP+2
          ENDDO
          IP=IP+NX
       ENDDO
            
       IP=NX+1
       DO J=1,NY2
          DO I=1,NX2
             OUT(I,J+NY/2)=WORK(IP)/FAC
             IP=IP+2
          ENDDO
          IP=IP+NX
       ENDDO
            
       IP=NX*NY+1
       DO J=1,NY2
          DO I=1,NX2
             OUT(I+NX/2,J)=WORK(IP)/FAC
             IP=IP+2
          ENDDO
          IP=IP+NX
       ENDDO
            
       IP=NX*NY+NX+1
       DO J=1,NY2
          DO I=1,NX2
             OUT(I,J)=WORK(IP)/FAC
             IP=IP+2
          ENDDO
          IP=IP+NX
       ENDDO
      ELSE
          
C 1D case
       DO I=1,NX2
         OUT(I+NX2,1)=WORK(IP)/FAC
         IP=IP+2
       ENDDO
            
       DO I=1,NX2
         OUT(I,1)=WORK(IP)/FAC
         IP=IP+2
       ENDDO
      ENDIF
           
      RETURN
      END
         
      SUBROUTINE DFOURT(DATA,NN,NDIM,ISIGN,IFORM,WORK)
C
C Double Precision version of FOURT.F
C
C Richard Hook, ST-ECF, October 1992
C
C-------------------Start of original code (modifed to be more readable)
C
C  June 30, 1983
C
C  Watch out for INTEGER*4 on VAX!! NN, ISIGN and IFORM must all
C  be dimensioned INTEGER*4 
C
C     The Fast Fourier Transform in USASI basic FORTRAN
C     -------------------------------------------------
C
C     TRANSFORM(J1,J2,...) = SUM(DATA(I1,I2,...)*W1**((I1-1)*(J1-1))
C                                 *W2**((I2-1)*(J2-1))*...),
C
C     where I1 and J1 run from 1 to NN(1) and W1=EXP(ISIGN*2*PI*
C     SQRT(-1)/NN(1)), etc.  there is no limit on the dimensionality
C     (number of subscripts) of the data array.  If an inverse
C     transform (ISIGN=+1) is performed upon an array of transformed
C     (ISIGN=-1) data, the original data will reappear,
C     multiplied by NN(1)*NN(2)*...  The array of input data may be
C     real or complex, at the programmers option, with a saving of
C     up to forty per cent in running time for real over complex.
C     (for fastest transform of real data, NN(1) should be even.)
C     the transform values are always complex, and are returned in the
C     original array of data, replacing the input data.  the length
C     of each dimension of the data array may be any integer.  the
C     program runs faster on composite integers than on primes, and is
C     particularly fast on numbers rich in factors of two.
C
C     Timing is in fact given by the following formula.  Let NTOT be the
C     total number of points (real or complex) in the data array, that
C     is, NTOT=NN(1)*NN(2)*...  decompose NTOT into its prime factors,
C     such as 2**K2 * 3**K3 * 5**K5 * ...  Let SUM2 be the sum of all
C     the factors of two in NTOT, that is, SUM2 = 2*K2.  Let SUMF be
C     the sum of all other factors of NTOT, that is, SUMF = 3*K3+5*K5+..
C     The time taken by a multidimensional transform on these NTOT data
C     point add time = six microseconds), T = 3000 + NTOT*(600+40*SUM2+
C     is T = T0 + NTOT*(T1+T2*SUM2+T3*SUMF). On the CDC 3300 (floating
C     175*SUMF) microseconds on complex data.
C
C     Implementation of the definition by summation will run in a time
C     proportional to NTOT*(NN(1)+NN(2)+...). For highly composite NTOT
C     the savings offered by this program can be dramatic. A one-dimen-
C     sional array 4000 in length will be transformed in 4000*(600+
C     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 seconds versus about 4000*
C     4000*175 = 2800 seconds for the straightforward technique.
C
C     [Note added by Richard Hook, October 1992:  as a rough guide  ]
C     [the time taken for a real and imaginary FFT of a 512*512     ]
C     [array is about 8s on a SPARCStation 2                        ]         
C
C     The Fast Fourier algorithm places two restrictions upon the
C     nature of the data beyond the usual restriction that
C     the data form one cycle of a periodic function.  They are--
C     1.  The number of input data and the number of transform values
C     must be the same.
C     2. Considering the data to be in the time domain,
C     they must be equi-spaced at intervals of DT.  Further, the trans-
C     form values, considered to be in frequency space, will be equi-
C     spaced from 0 to 2*PI*(NN(I)-1)/(NN(I)*DT) at intervals of
C     2*PI/(NN(I)*DT) for each dimension of length NN(I).  Of course,
C     dt need not be the same for every dimension.
C
C     The calling sequence is--
C
C     CALL DFOURT(DATA,NN,NDIM,ISIGN,IFORM,WORK)
C
C     DATA is the array used to hold the real and imaginary parts
C     of the data on input and the transform values on output.  It
C     is a multidimensional floating point array, with the real and
C     imaginary parts of a datum stored immediately adjacent in storage
C     (such as FORTRAN IV places them).  The extent of each dimension
C     is given in the integer array NN, of length NDIM.  ISIGN is -1
C     to indicate a forward transform (exponential sign is -) and +1
C     for an inverse transform (sign is +).  IFORM is +1 if the data and
C     the transform values are complex.  It is 0 if the data are real
C     but the transform values are complex.  If it is 0, the imaginary
C     parts of the data should be set to zero.  As explained above, the
C     transform values are always complex and are stored in array data.
C     work is an array used for working storage.  It is not necessary
C     if all the dimensions of the data are powers of two.  In this case
C     it may be replaced by 0 in the calling sequence.  Thus, use of
C     powers of two can free a good deal of storage.  If any dimension
C     is not a power of two, this array must be supplied.  It is
C     floating point, one dimensional of length equal to twice the
C     largest array dimension (i.e., NN(I) ) that is not a power of
C     two.  Therefore, in one dimension for a non power of two,
C     work occupies as many storage locations as data.  If supplied,
C     work must not be the same array as data.  All subscripts of all
C     arrays begin at one.
C
C     Example 1.  Three-dimensional forward Fourier Transform of a
C     complex array dimensioned 32 by 25 by 13 in FORTRAN IV.
C
C     DOUBLE PRECISION DATA(32,25,13),WORK(50),NN(3)
C     COMPLEX DATA
C     DATA NN/32,25,13/
C     DO 1 I=1,32
C     DO 1 J=1,25
C     DO 1 K=1,13
C  1  DATA(I,J,K)=COMPLEX VALUE
C     CALL DFOURT(DATA,NN,3,-1,1,WORK)
C
C     Example 2.  One-dimensional forward transform of a real array of
C     length 64 in FORTRAN II.
C
C     DOUBLE PRECISION DATA(2,64)
C     DO 2 I=1,64
C     DATA(1,I)=REAL PART
C  2  DATA(2,I)=0.
C     CALL DFOURT(DATA,64,1,-1,0,0)
C
C     There are no error messages or error halts in this program.  The
C     program returns immediately if ndim or any nn(i) is less than one.
C
C     Program by Norman Brenner from the BASIC program by Charles
C     Rader (both of MIT Lincoln Laboratory).  May 1967.  The idea
C     for the digit reversal was suggested by Ralph Alter (also MIT LL).
C     This is the fastest and most versatile version of the FFT known
C     to the author.  A program called four2 is available that also
C     performs the fast fourier transform and is written in usasi basic
C     fortran.  It is about one third as long and restricts the
C     dimensions of the input array (which must be complex) to be powers
C     of two.  Another program, called four1, is one tenth as long and
C     runs two thirds as fast on a one-dimensional complex array whose
C     length is a power of two.
C
C     Reference--
C     Fast Fourier Transforms for Fun and Profit, W. Gentleman and
C     G. Sande, 1966 Fall Joint Computer Conference.
C
C     The work reported in this document was performed at Lincoln Lab-
C     oratory, a center for research operated by Massachusetts Institute
C     of Technology, with the support of the U.S. Air Force under
C     contract af 19(628)-5167.
C
C-------------End of original comment block.
C
      IMPLICIT  NONE
C
      DOUBLE PRECISION DATA(*),WORK(*)
      DOUBLE PRECISION TWOPI,RTHLF,THETA,WR,WI
      DOUBLE PRECISION W2R,W2I,W3R,W3I,SUMR,SUMI,DIFR,DIFI
      DOUBLE PRECISION TEMPR,TEMPI,T2R,T2I,T3R,T3I,WTEMP
      DOUBLE PRECISION T4R,T4I,WSTPR,WSTPI,THETM,WMINR,WMINI
      DOUBLE PRECISION TWOWR,SR,SI,OLDSR,OLDSI,STMPR,STMPI
      DOUBLE PRECISION U1R,U1I,U2R,U2I,U3R,U3I,U4I,U4R

      INTEGER IFACT(32),NN(1)
      INTEGER ISIGN,IFORM
      INTEGER NDIM,NTOT,NHALF,IMAX,IMIN,JMAX
      INTEGER NP1,NP2,I,J,N,M,NTWO,IF,IDIV,IQUOT,IREM,INON2,NON2P
      INTEGER ICASE,IFMIN,I1RNG,IDIM,NP0,NPREV
      INTEGER I1,I2,I3,J1,J2,J3,K1,K2,K3,IPAR,KDIF,NP1TW,MMAX
      INTEGER NP2HF,I1MAX,NWORK,IFP1,IFP2,I2MAX,LMAX,L
      INTEGER KMIN,KSTEP,K4,J2MAX,JMIN,J3MAX

C Initialise constants
      TWOPI=6.28318530717959D0
      RTHLF=DSQRT(2.0D0)/2.0D0

C Start of executable code
      IF(NDIM-1)920,1,1
1     NTOT=2
      DO 2 IDIM=1,NDIM
      IF(NN(IDIM))920,920,2
2     NTOT=NTOT*NN(IDIM)
C
C     main loop for each dimension
C
      NP1=2
      DO 910 IDIM=1,NDIM
      N=NN(IDIM)
      NP2=NP1*N
      IF(N-1)920,900,5
C
C     is n a power of two and if not, what are its factors
C
5     M=N
      NTWO=NP1
      IF=1
      IDIV=2
10    IQUOT=M/IDIV
      IREM=M-IDIV*IQUOT
      IF(IQUOT-IDIV)50,11,11
11    IF(IREM)20,12,20
12    NTWO=NTWO+NTWO
      IFACT(IF)=IDIV
      IF=IF+1
      M=IQUOT
      GO TO 10
20    IDIV=3
      INON2=IF
30    IQUOT=M/IDIV
      IREM=M-IDIV*IQUOT
      IF(IQUOT-IDIV)60,31,31
31    IF(IREM)40,32,40
32    IFACT(IF)=IDIV
      IF=IF+1
      M=IQUOT
      GO TO 30
40    IDIV=IDIV+2
      GO TO 30
50    INON2=IF
      IF(IREM)60,51,60
51    NTWO=NTWO+NTWO
      GO TO 70
60    IFACT(IF)=M
70    NON2P=NP2/NTWO
C
C     Separate four cases--
C        1. Complex transform
C        2. Real transform for the 2nd, 3rd, etc. dimension.  method--
C           transform half the data, supplying the other half by con-
C           jugate symmetry.
C        3. Real transform for the 1st dimension, n odd.  method--
C           set the imaginary parts to zero.
C        4. Real transform for the 1st dimension, n even.  method--
C           transform a complex array of length n/2 whose real parts
C           are the even numbered real values and whose imaginary parts
C           are the odd numbered real values.  separate and supply
C           the second half by conjugate symmetry.
C
      ICASE=1
      IFMIN=1
      I1RNG=NP1
      IF(IDIM-4)74,100,100
74    IF(IFORM)71,71,100
71    ICASE=2
      I1RNG=NP0*(1+NPREV/2)
      IF(IDIM-1)72,72,100
72    ICASE=3
      I1RNG=NP1
      IF(NTWO-NP1)100,100,73
73    ICASE=4
      IFMIN=2
      NTWO=NTWO/2
      N=N/2
      NP2=NP2/2
      NTOT=NTOT/2
      I=1
      DO 80 J=1,NTOT
      DATA(J)=DATA(I)
80    I=I+2
C
C     Shuffle data by bit reversal, since n=2**k.  as the shuffling
C     can be done by simple interchange, no working array is needed
C
100   IF(NON2P-1)101,101,200
101   NP2HF=NP2/2
      J=1
      DO 150 I2=1,NP2,NP1
      IF(J-I2)121,130,130
121   I1MAX=I2+NP1-2
      DO 125 I1=I2,I1MAX,2
      DO 125 I3=I1,NTOT,NP2
      J3=J+I3-I2
      TEMPR=DATA(I3)
      TEMPI=DATA(I3+1)
      DATA(I3)=DATA(J3)
      DATA(I3+1)=DATA(J3+1)
      DATA(J3)=TEMPR
125   DATA(J3+1)=TEMPI
130   M=NP2HF
140   IF(J-M)150,150,141
141   J=J-M
      M=M/2
      IF(M-NP1)150,140,140
150   J=J+M
      GO TO 300
C
C     Shuffle data by digit reversal for general n
C
200   NWORK=2*N
      DO 270 I1=1,NP1,2
      DO 270 I3=I1,NTOT,NP2
      J=I3
      DO 260 I=1,NWORK,2
      IF(ICASE-3)210,220,210
210   WORK(I)=DATA(J)
      WORK(I+1)=DATA(J+1)
      GO TO 240
220   WORK(I)=DATA(J)
      WORK(I+1)=0.0D0
240   IFP2=NP2
      IF=IFMIN
250   IFP1=IFP2/IFACT(IF)
      J=J+IFP1
      IF(J-I3-IFP2)260,255,255
255   J=J-IFP2
      IFP2=IFP1
      IF=IF+1
      IF(IFP2-NP1)260,260,250
260   CONTINUE
      I2MAX=I3+NP2-NP1
      I=1
      DO 270 I2=I3,I2MAX,NP1
      DATA(I2)=WORK(I)
      DATA(I2+1)=WORK(I+1)
270   I=I+2
C
C     Main loop for factors of two.
C     W=EXP(ISIGN*2*PI*SQRT(-1)*M/(4*MMAX)).  Check for W=ISIGN*SQRT(-1)
C     and repeat for W=W*(1+ISIGN*SQRT(-1))/SQRT(2).
C
300   IF(NTWO-NP1)600,600,305
305   NP1TW=NP1+NP1
      IPAR=NTWO/NP1
310   IF(IPAR-2)350,330,320
320   IPAR=IPAR/4
      GO TO 310
330   DO 340 I1=1,I1RNG,2
      DO 340 K1=I1,NTOT,NP1TW
      K2=K1+NP1
      TEMPR=DATA(K2)
      TEMPI=DATA(K2+1)
      DATA(K2)=DATA(K1)-TEMPR
      DATA(K2+1)=DATA(K1+1)-TEMPI
      DATA(K1)=DATA(K1)+TEMPR
340   DATA(K1+1)=DATA(K1+1)+TEMPI
350   MMAX=NP1
360   IF(MMAX-NTWO/2)370,600,600
370   LMAX=MAX0(NP1TW,MMAX/2)
      DO 570 L=NP1,LMAX,NP1TW
      M=L
      IF(MMAX-NP1)420,420,380
380   THETA=-TWOPI*DFLOAT(L)/DFLOAT(4*MMAX)
      IF(ISIGN)400,390,390
390   THETA=-THETA
400   WR=DCOS(THETA)
      WI=DSIN(THETA)
410   W2R=WR*WR-WI*WI
      W2I=2.0D0*WR*WI
      W3R=W2R*WR-W2I*WI
      W3I=W2R*WI+W2I*WR
420   DO 530 I1=1,I1RNG,2
      KMIN=I1+IPAR*M
      IF(MMAX-NP1)430,430,440
430   KMIN=I1
440   KDIF=IPAR*MMAX
450   KSTEP=4*KDIF
      IF(KSTEP-NTWO)460,460,530
460   DO 520 K1=KMIN,NTOT,KSTEP
      K2=K1+KDIF
      K3=K2+KDIF
      K4=K3+KDIF
      IF(MMAX-NP1)470,470,480
470   U1R=DATA(K1)+DATA(K2)
      U1I=DATA(K1+1)+DATA(K2+1)
      U2R=DATA(K3)+DATA(K4)
      U2I=DATA(K3+1)+DATA(K4+1)
      U3R=DATA(K1)-DATA(K2)
      U3I=DATA(K1+1)-DATA(K2+1)
      IF(ISIGN)471,472,472
471   U4R=DATA(K3+1)-DATA(K4+1)
      U4I=DATA(K4)-DATA(K3)
      GO TO 510
472   U4R=DATA(K4+1)-DATA(K3+1)
      U4I=DATA(K3)-DATA(K4)
      GO TO 510
480   T2R=W2R*DATA(K2)-W2I*DATA(K2+1)
      T2I=W2R*DATA(K2+1)+W2I*DATA(K2)
      T3R=WR*DATA(K3)-WI*DATA(K3+1)
      T3I=WR*DATA(K3+1)+WI*DATA(K3)
      T4R=W3R*DATA(K4)-W3I*DATA(K4+1)
      T4I=W3R*DATA(K4+1)+W3I*DATA(K4)
      U1R=DATA(K1)+T2R
      U1I=DATA(K1+1)+T2I
      U2R=T3R+T4R
      U2I=T3I+T4I
      U3R=DATA(K1)-T2R
      U3I=DATA(K1+1)-T2I
      IF(ISIGN)490,500,500
490   U4R=T3I-T4I
      U4I=T4R-T3R
      GO TO 510
500   U4R=T4I-T3I
      U4I=T3R-T4R
510   DATA(K1)=U1R+U2R
      DATA(K1+1)=U1I+U2I
      DATA(K2)=U3R+U4R
      DATA(K2+1)=U3I+U4I
      DATA(K3)=U1R-U2R
      DATA(K3+1)=U1I-U2I
      DATA(K4)=U3R-U4R
520   DATA(K4+1)=U3I-U4I
      KDIF=KSTEP
      KMIN=4*(KMIN-I1)+I1
      GO TO 450
530   CONTINUE
      M=M+LMAX
      IF(M-MMAX)540,540,570
540   IF(ISIGN)550,560,560
550   TEMPR=WR
      WR=(WR+WI)*RTHLF
      WI=(WI-TEMPR)*RTHLF
      GO TO 410
560   TEMPR=WR
      WR=(WR-WI)*RTHLF
      WI=(TEMPR+WI)*RTHLF
      GO TO 410
570   CONTINUE
      IPAR=3-IPAR
      MMAX=MMAX+MMAX
      GO TO 360
C
C     Main loop for factors not equal to two.
C     W=EXP(ISIGN*2*PI*SQRT(-1)*(J1+J2-I3-1)/IFP2)
C
600   IF(NON2P-1)700,700,601
601   IFP1=NTWO
      IF=INON2
610   IFP2=IFACT(IF)*IFP1
      THETA=-TWOPI/DFLOAT(IFACT(IF))
      IF(ISIGN)612,611,611
611   THETA=-THETA
612   WSTPR=DCOS(THETA)
      WSTPI=DSIN(THETA)
      DO 650 J1=1,IFP1,NP1
      THETM=-TWOPI*DFLOAT(J1-1)/DFLOAT(IFP2)
      IF(ISIGN)614,613,613
613   THETM=-THETM
614   WMINR=DCOS(THETM)
      WMINI=DSIN(THETM)
      I1MAX=J1+I1RNG-2
      DO 650 I1=J1,I1MAX,2
      DO 650 I3=I1,NTOT,NP2
      I=1
      WR=WMINR
      WI=WMINI
      J2MAX=I3+IFP2-IFP1
      DO 640 J2=I3,J2MAX,IFP1
      TWOWR=WR+WR
      J3MAX=J2+NP2-IFP2
      DO 630 J3=J2,J3MAX,IFP2
      JMIN=J3-J2+I3
      J=JMIN+IFP2-IFP1
      SR=DATA(J)
      SI=DATA(J+1)
      OLDSR=0.0D0
      OLDSI=0.0D0
      J=J-IFP1
620   STMPR=SR
      STMPI=SI
      SR=TWOWR*SR-OLDSR+DATA(J)
      SI=TWOWR*SI-OLDSI+DATA(J+1)
      OLDSR=STMPR
      OLDSI=STMPI
      J=J-IFP1
      IF(J-JMIN)621,621,620
621   WORK(I)=WR*SR-WI*SI-OLDSR+DATA(J)
      WORK(I+1)=WI*SR+WR*SI-OLDSI+DATA(J+1)
630   I=I+2
      WTEMP=WR*WSTPI
      WR=WR*WSTPR-WI*WSTPI
640   WI=WI*WSTPR+WTEMP
      I=1
      DO 650 J2=I3,J2MAX,IFP1
      J3MAX=J2+NP2-IFP2
      DO 650 J3=J2,J3MAX,IFP2
      DATA(J3)=WORK(I)
      DATA(J3+1)=WORK(I+1)
650   I=I+2
      IF=IF+1
      IFP1=IFP2
      IF(IFP1-NP2)610,700,700
C
C     Complete a real transform in the 1st dimension, N even, by con-
C     jugate symmetries.
C
700   GO TO (900,800,900,701),ICASE
701   NHALF=N
      N=N+N
      THETA=-TWOPI/DFLOAT(N)
      IF(ISIGN)703,702,702
702   THETA=-THETA
703   WSTPR=DCOS(THETA)
      WSTPI=DSIN(THETA)
      WR=WSTPR
      WI=WSTPI
      IMIN=3
      JMIN=2*NHALF-1
      GO TO 725
710   J=JMIN
      DO 720 I=IMIN,NTOT,NP2
      SUMR=(DATA(I)+DATA(J))/2.0D0
      SUMI=(DATA(I+1)+DATA(J+1))/2.0D0
      DIFR=(DATA(I)-DATA(J))/2.0D0
      DIFI=(DATA(I+1)-DATA(J+1))/2.0D0
      TEMPR=WR*SUMI+WI*DIFR
      TEMPI=WI*SUMI-WR*DIFR
      DATA(I)=SUMR+TEMPR
      DATA(I+1)=DIFI+TEMPI
      DATA(J)=SUMR-TEMPR
      DATA(J+1)=-DIFI+TEMPI
720   J=J+NP2
      IMIN=IMIN+2
      JMIN=JMIN-2
      WTEMP=WR*WSTPI
      WR=WR*WSTPR-WI*WSTPI
      WI=WI*WSTPR+WTEMP
725   IF(IMIN-JMIN)710,730,740
730   IF(ISIGN)731,740,740
731   DO 735 I=IMIN,NTOT,NP2
735   DATA(I+1)=-DATA(I+1)
740   NP2=NP2+NP2
      NTOT=NTOT+NTOT
      J=NTOT+1
      IMAX=NTOT/2+1
745   IMIN=IMAX-2*NHALF
      I=IMIN
      GO TO 755
750   DATA(J)=DATA(I)
      DATA(J+1)=-DATA(I+1)
755   I=I+2
      J=J-2
      IF(I-IMAX)750,760,760
760   DATA(J)=DATA(IMIN)-DATA(IMIN+1)
      DATA(J+1)=0.
      IF(I-J)770,780,780
765   DATA(J)=DATA(I)
      DATA(J+1)=DATA(I+1)
770   I=I-2
      J=J-2
      IF(I-IMIN)775,775,765
775   DATA(J)=DATA(IMIN)+DATA(IMIN+1)
      DATA(J+1)=0.
      IMAX=IMIN
      GO TO 745
780   DATA(1)=DATA(1)+DATA(2)
      DATA(2)=0.
      GO TO 900
C
C     Complete a real transform for the 2nd, 3rd, etc. dimension by
C     conjugate symmetries.
C
800   IF(I1RNG-NP1)805,900,900
805   DO 860 I3=1,NTOT,NP2
      I2MAX=I3+NP2-NP1
      DO 860 I2=I3,I2MAX,NP1
      IMAX=I2+NP1-2
      IMIN=I2+I1RNG
      JMAX=2*I3+NP1-IMIN
      IF(I2-I3)820,820,810
810   JMAX=JMAX+NP2
820   IF(IDIM-2)850,850,830
830   J=JMAX+NP0
      DO 840 I=IMIN,IMAX,2
      DATA(I)=DATA(J)
      DATA(I+1)=-DATA(J+1)
840   J=J-2
850   J=JMAX
      DO 860 I=IMIN,IMAX,NP0
      DATA(I)=DATA(J)
      DATA(I+1)=-DATA(J+1)
860   J=J-NP0
C
C     End of loop on each dimension
C
900   NP0=NP1
      NP1=NP2
910   NPREV=N
920   RETURN
      END

      SUBROUTINE DIRCON(DATA,NX,NY,PSF,PNX,PNY,OUT,FLAG)
C
C DIRCON - direct convolution, not using FFTs.
C
C First simple version, Richard Hook, June 1993
C
      IMPLICIT NONE

      INTEGER NX,NY,PNX,PNY,FLAG
      DOUBLE PRECISION DATA(NX,NY),PSF(PNX,PNY),OUT(NX,NY)

C Local variables
      INTEGER I,J,N,M,IS,JS
      DOUBLE PRECISION T

      DO J=PNY/2,NY-PNY/2
         JS=J-PNY/2
         DO I=PNX/2,NX-PNX/2
            IS=I-PNX/2

            T=0.0D0
            DO N=1,PNY
               DO M=1,PNX
                  T=T+DATA(IS+M,JS+N)*PSF(M,N)
               ENDDO
            ENDDO

            OUT(I,J)=T
         ENDDO
      ENDDO

      END

      SUBROUTINE SIRCON(DATA,NX,NY,PSF,PNX,PNY,OUT,FLAG)
C
C SIRCON - direct convolution, not using FFTs.
C Single precision version.
C
C First simple version, Richard Hook, June 1993
C
      IMPLICIT NONE

      INTEGER NX,NY,PNX,PNY,FLAG
      REAL DATA(NX,NY),PSF(PNX,PNY),OUT(NX,NY)

C Local variables
      INTEGER I,J,N,M,IS,JS
      REAL T

      DO J=PNY/2,NY-PNY/2
         JS=J-PNY/2
         DO I=PNX/2,NY-PNX/2
            IS=I-PNX/2

            T=0.0
            DO N=1,PNY
               DO M=1,PNX
                  T=T+DATA(IS+M,JS+N)*PSF(M,N)
               ENDDO
            ENDDO

            OUT(I,J)=T
         ENDDO
      ENDDO

      END

      SUBROUTINE FILCON(DATA,NX,NY,VAL)
C
C Fill up a double precision array with a constant value
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      DOUBLE PRECISION DATA(NX,NY),VAL

      DO J=1,NY
         DO I=1,NX
            DATA(I,J)=VAL
          ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SILCON(DATA,NX,NY,VAL)
C
C Fill up a single precision array with a constant value
C
      IMPLICIT NONE

      INTEGER NX,NY,I,J
      REAL DATA(NX,NY),VAL

      DO J=1,NY
         DO I=1,NX
            DATA(I,J)=VAL
          ENDDO
      ENDDO

      RETURN
      END


