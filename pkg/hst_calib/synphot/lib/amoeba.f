C*AMOEBA ... downhill simplex fitting program
C+
	SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,FUNK,ITER)
*
* adjust parameters to minimize value of a function
*
* input:
*	P	R4(MP,NP)	simplex vertices
*	Y	R4(NDIM+1)	function value at each vertex
*	MP	I4		max number of vertices
*	NP	I4		max number of parameters
*	NDIM	I4		actual number of parameters
*	FTOL	R4		relative tolerance for convergence
*	ITER	I4		maximum iterations
* output:
*	P	R4(MP,NP)	final vertices of simplex
*	Y	R4(NDIM+1)	FUNK values at each vertex
*	ITER	I4		number of iterations
C--
*
* May 1989 Keith Horne @ STScI - adapted from Numerical Recipes
* May 1989 KDH @ STScI - modified convergence test ala Kip Kuntz.
* Jun 1989 KDH @ STScI - convergence messages to terminal.
*
	PARAMETER (NMAX=100)
	DIMENSION P(MP,NP),Y(MP),PR(NMAX),PRR(NMAX),PBAR(NMAX)
	LOGICAL CONTRACT
	EXTERNAL FUNK
C
	ALPHA=1.0
	BETA=0.5
	GAMMA=2.0
	ITMAX = ITER
	MPTS=NDIM+1
	ITER=0
1	ILO=1
	IF(Y(1).GT.Y(2))THEN
		IHI=1
		INHI=2
	ELSE
		IHI=2
		INHI=1
	ENDIF
	DO I=1,MPTS
		IF(Y(I).LT.Y(ILO)) ILO=I
		IF(Y(I).GT.Y(IHI))THEN
			INHI=IHI
			IHI=I
		ELSE IF(Y(I).GT.Y(INHI))THEN
			IF(I.NE.IHI) INHI=I
		ENDIF
	ENDDO

	RTOL=2.*ABS(Y(IHI)-Y(ILO))/(ABS(Y(IHI))+ABS(Y(ILO)))

* converge if chi^2 range is small last operation is contraction

      IF( RTOL.LT.FTOL ) then
      IF( CONTRACT ) THEN
	WRITE(*,*) 'Amoeba converged. Function range', RTOL, ' <', FTOL
	GOTO 1000
      END IF

* converge if precision limit reached in all parameters

	TMAX = 0
	PRECISE = 1E-7
      DO I=1,NDIM
	TEST = ABS( P(I,IHI) ) + ABS( P(I,ILO) )
      IF( TEST.GE. PRECISE ) THEN
	TEST = ABS( P(I,IHI) - P(I,ILO) )
     #		/ ( ABS( P(I,IHI) ) + ABS( P(I,ILO) ) )
      END IF	
	TMAX = MAX( TMAX, TEST )
      END DO
      IF( TMAX .LT. PRECISE ) THEN
	WRITE(*,*) 'Amoeba converged. Parameter range', TMAX, ' <', PRECISE
	GOTO 1000
      END IF
      END IF

* converge if maximum iterations exceeded

      IF(ITER.EQ.ITMAX) THEN
      IF( ITMAX.GT.1 ) THEN
	WRITE(*,*) 'Amoeba iteration count exceeded.', ITER
      END IF
	GOTO 1000
      END IF
	CONTRACT = .FALSE.
	ITER=ITER+1
	DO J=1,NDIM
		PBAR(J)=0.
	ENDDO
	DO I=1,MPTS
		IF(I.NE.IHI)THEN
			DO J=1,NDIM
				PBAR(J)=PBAR(J)+P(I,J)
			ENDDO
		ENDIF
	ENDDO
	DO J=1,NDIM
		PBAR(J)=PBAR(J)/NDIM
		PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
	ENDDO
	YPR=FUNK(PR)
	IF(YPR.LE.Y(ILO))THEN
		DO J=1,NDIM
			PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
		ENDDO
		YPRR=FUNK(PRR)
		IF(YPRR.LT.Y(ILO))THEN
			DO J=1,NDIM
				P(IHI,J)=PRR(J)
			ENDDO
			Y(IHI)=YPRR
		ELSE
			DO J=1,NDIM
				P(IHI,J)=PR(J)
			ENDDO
			Y(IHI)=YPR
		ENDIF
	ELSE IF(YPR.GE.Y(INHI))THEN
		IF(YPR.LT.Y(IHI))THEN
			DO J=1,NDIM
				P(IHI,J)=PR(J)
			ENDDO
			Y(IHI)=YPR
		ENDIF
		DO J=1,NDIM
			PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
		ENDDO
		YPRR=FUNK(PRR)
		IF(YPRR.LT.Y(IHI))THEN
			DO J=1,NDIM
				P(IHI,J)=PRR(J)
			ENDDO
			Y(IHI)=YPRR
		ELSE
			contract = .true.
			DO I=1,MPTS
				IF(I.NE.ILO)THEN
					DO J=1,NDIM
						PR(J)=0.5*(P(I,J)+P(ILO,J))
						P(I,J)=PR(J)
					ENDDO
					Y(I)=FUNK(PR)
				ENDIF
			ENDDO
		ENDIF
	ELSE
		DO J=1,NDIM
			P(IHI,J)=PR(J)
		ENDDO
		Y(IHI)=YPR
	ENDIF
	GO TO 1

* normal return

 1000	RETURN
	END
