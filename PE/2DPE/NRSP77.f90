MODULE NRSP77
 	USE NRTYPE
	USE NRUTIL
 	IMPLICIT NONE
 	!-----------------------------------------------------------
 	CONTAINS
	!-----------------------------------------------------------
 	SUBROUTINE DSPRSIN(A,N,Np,THRESH,Nmax,SA,IJA)
  	INTEGER(I4B) :: N, Nmax, Np, IJA(Nmax)
  	INTEGER(I4B) :: I, J, K
  	REAL(DP) :: THRESH
  	COMPLEX(DP) :: A(Np,Np), SA(Nmax)
  	DO J = 1,N
   	SA(J) = A(J,J)
  	END DO
  	IJA(1) = N+2
  	K = N+1
 	  DO I = 1,N
   		DO J = 1,N
    			IF (ABS(A(I,J)).GT.THRESH) THEN
    				IF (I.NE.J) THEN
     					K = K+1
      					IF (K.GT.Nmax) THEN
       						WRITE(*,*) 'Nmax too small. Press Enter to continue'
      						READ(*,*)
      					END IF
      					SA(K) = A(I,J)
      					IJA(K) = j
     				END IF
    			END IF
   		END DO
   		IJA(I+1) = K+1
  	END DO
 	END SUBROUTINE 
 	!-----------------------------------------------------------
 	SUBROUTINE DSPRSAX(SA,IJA,X,B,N)
 	INTEGER(I4B) :: N, IJA(:)
  	INTEGER(I4B) :: I, K
  	COMPLEX(DP) :: B(N), X(N), SA(:)
  	IF (IJA(1).NE.N+2) THEN
   		WRITE( *, * ) 'Mismatched vector and matrix. Press Enter to continue'
   		READ( *, * )
  	END IF
  	DO I = 1,N
   		B(I) = SA(I)*X(I)
   		DO K = IJA(I),IJA(I+1)-1
    			B(I) = B(I)+SA(K)*X(IJA(K))
   		END DO
  	END DO
 	END SUBROUTINE 
 	!-----------------------------------------------------------
 	SUBROUTINE DSPRSTX(SA,IJA,X,B,N)
  	INTEGER(I4B) :: N, IJA(:)
  	INTEGER(I4B) :: I, J, K
  	COMPLEX(DP) :: B(N), X(N), SA(:)
  	IF (IJA(1).NE.N+2) THEN
   		WRITE( *, * ) 'Mismatched vector and matrix. Press Enter to continue'
   		READ( *, * )
  	END IF
  	DO I = 1,N
   		B(I) = SA(I)*X(I)
   		DO K = IJA(I),IJA(I+1)-1
    			J = IJA(K)
    			B(J) = B(J)+SA(K)*X(IJA(K))
   		END DO
  	END DO
 	END SUBROUTINE 
 	!-----------------------------------------------------------
 	SUBROUTINE DSPRSTP(SA,IJA,SB,IJB)
  	INTEGER(I4B) :: IJA(:), IJB(:)
  	COMPLEX(DP) :: SA(:), SB(:), V
  	INTEGER(I4B) :: J, JL, JM, JP, JU, K, M, N2, NOFF, INC, IV
  	N2 = IJA(1)
  	DO J = 1,N2-2
   		SB(J) = SA(J)
  	END DO
  	CALL IINDEXX(IJA(N2-1)-IJA(1),IJA(N2),IJB(N2))
  	JP = 0
  	DO K = IJA(1),IJA(N2-1)-1
   		M = IJB(K)+N2-1
   		SB(K) = SA(M)
   		DO J = JP+1,IJA(M)
    		IJB(J) = K
   		END DO
   		JP = IJA(M)
   		JL = 1
   		JU = N2-1
   		DO WHILE (JU-JL.GT.1) THEN
    			JM = (JU+JL)/2
    			IF (IJA(JM).GT.M) THEN
     				JU = JM
    			ELSE
     				JL = JM
    			END IF
   		END DO
   		IJB(K) = JL
  	END DO
  	DO J = JP+1,N2-1
  		IJB(J) = IJA(N2-1)
  	END DO
  	DO J = 1,N2-2
   		JL = IJB(J+1)-IJB(J)
   		NOFF = IJB(J)-1
   		INC = 1
1  		INC = 3*INC+1
   		IF (INC.LE.JL) GOTO 1
   		WHILE DO (INC.GT.1)
    			INC = INC/3
    			DO K = NOFF+INC+1,NOFF+JL

   	END DO 
 	END SUBROUTINE
 	!-----------------------------------------------------------
 	SUBROUTINE DSPRSTM(SA,IJA,SB,IJB,THRESH,NMAX,SC,IJC)
  	INTEGER(I4B) :: NMAX, IJA(:), IJB(:), IJC(NMAX)
  	REAL(DP) :: THRESH
  	COMPLEX(DP) :: SA(:), SB(:), SC(NMAX), SUMM
  	INTEGER(I4B) :: I, IJMA, IJMB, J, K, MA, MB, MBB
  	IF (IJA(1).NE.IJB(1)) THEN
   		WRITE(*,*) 'Size does not match. Press Enter to continue'
   		READ(*,*)
  	END IF
  	K = IJA(1)
  	IJC(1) = K
  	DO I = 1,IJA(1)-2
  		DO J = 1,IJB(1)-2
    			IF (I.EQ.J) THEN
     				SUMM = SA(I)*SB(J)
    			ELSE
     				SUMM = 0.d0
    			END IF
    			MB = IJB(J)
    			DO MA = IJA(I),IJA(I+1)-1
     				IJMA = IJA(MA)
     				IF (IJMA.EQ.J) THEN
      					SUMM = SUMM+SA(MA)*SB(J)
     				ELSE
      					DO WHILE (MB.LT.IJB(J+1))
       						IJMB = IJB(MB)
       						IF (IJMB.EQ.I) THEN
        						SUMM = SUMM+SA(I)*SB(MB)
        						MB = MB+1
       						ELSEIF (IJMB.LT.IJMA) THEN
        						MB = MB+1
       						ELSEIF (IJMB.EQ.IJMA) THEN
        						SUMM = SUMM+SA(MA)*SB(MB)
        						MB = MB+1
       						END IF
      					END DO
     				END IF
    			END DO
    			DO MBB = MB,IJB(J+1)-1
     				IF (IJB(MBB).EQ.I) THEN
      					SUMM = SUMM + SA(I)*SB(MBB)
     				END IF
    			END DO
    			IF (I.EQ.J) THEN
     				SC(I) = SUMM
    			ELSE IF (ABS(SUMM).GT.THRESH) THEN
     				IF (K.GT. NMAX) THEN
      					WRITE(*,*) 'DSPRSTM : NMAX too small. Press Enter to continue'
      					READ(*,*)
     				END IF
     				SC(K) = SUMM
     				IJC(K) = J
     				K = K+1
    			END IF
   		END DO
   		IJC(I+1) = K
  	END DO
 	END SUBROUTINE 
 	!-----------------------------------------------------------
END MODULE NRSP77