      SUBROUTINE MEDFIL(IN,OUT,NX,NY)
C
C A simple median filter using a 3 by 3 box
C
      IMPLICIT NONE

      INTEGER NX,NY
      DOUBLE PRECISION IN(NX,NY),OUT(NX,NY),SUB(9)

      INTEGER I,J,M,N,NN

      DO J=2,NY-1
         DO I=2,NX-1
            NN=1
            DO N=-1,1
               DO M=-1,1
                  SUB(NN)=IN(I+M,J+N)
                  NN=NN+1
               ENDDO
            ENDDO
 
            CALL SORT(SUB,9)

            OUT(I,J)=SUB(5)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE SORT(RA,N)
      INTEGER N,I,J,L,IR
      DOUBLE PRECISION RA(N),RRA
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
