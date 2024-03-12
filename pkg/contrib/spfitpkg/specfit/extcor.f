        REAL FUNCTION EXTCOR(WAVE,EBV)
	REAL WAVE, EBV
C
C      MEAN GALACTIC CURVE (Seaton,M. 1979, MNRAS,187, 73p)
C      THE METHOD EMPLOYED IS LINEAR INTERPOLATION.

       REAL MAG1, M, X(18), MAG(18)
       REAL DY, DX, B, Z

C      Data for mean galactic with x less than 2.7
       DATA X /1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,
     &  2.3,2.4,2.5,2.6,2.7/
       DATA MAG /1.36,1.64,1.84,2.04,2.24,2.44,2.66,2.88,3.14,3.36,
     &  3.56,3.77,3.96,4.15,4.26,4.40,4.52,4.64/

       Z = 10000. / WAVE

C      GAK modification for a Herschel 36 w/o bump from Fitzpatrick and Massa
C      1988, ApJ, 328, 734.  Weak bump relative to Seaton, slightly stronger
C      FUV upturn.
C	mag1 = 3.14 + 1.680 + 0.055 * Z
C	if ( z .ge. 5.9 )
C     &    mag1=mag1+0.312*(0.5392*(Z-5.9)**2+0.0564*(z-5.9)**3)
C       extcor = 10. ** ( -0.4 * EBV * MAG1)
C       return

C	GAK modification for Theta1 Ori C from Fitspatrick and Mass (1988).
C	This star has weak bump and weak UV upturn relative to Seaton curve.
C	drude = z**2 / ( (z**2 - 4.635**2)**2 + (0.846 * z)**2)
C	mag1 = 3.14 + 1.251 + 0.033 * Z + 1.331 * drude
C	if ( z .ge. 5.9 )
C     &    mag1=mag1+0.186*(0.5392*(Z-5.9)**2+0.0564*(z-5.9)**3)
C       extcor = 10. ** ( -0.4 * EBV * MAG1)
C       return

            MAG1 = 0.0
            IF( Z .LE. 1.0 ) GO TO 90

            IF( Z .GE. 2.7 ) GO TO 30

              DO 10 I = 1,18
   10               IF (X(I).GT. Z) GOTO 20
   20         DY = MAG(I) - MAG(I - 1)
              DX = X(I) - X(I - 1)
              M = DY / DX 
              B = MAG(I) - M * X(I)
              MAG1 = M * Z + B
            GO TO 90

   30             IF( Z .GE. 3.65 ) GO TO 50
            MAG1=1.56+1.048*Z+1.01/((Z-4.60)**2+0.280)
            GO TO 90

   50             IF( Z .GE. 7.14 ) GO TO 70
            MAG1=2.29+0.848*Z+1.01/((Z-4.60)**2+0.280)
            GO TO 90

   70             MAG1=16.17-3.20*Z+0.2975*Z**2
              GO TO 90

   90             EXTCOR = 10. ** ( -.4 * MAG1 * EBV )

       RETURN
       END
