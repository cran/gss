      SUBROUTINE DCHDC(A,LDA,P,WORK,JPVT,JOB,INFO)
      INTEGER LDA,P,JPVT(1),JOB,INFO
      DOUBLE PRECISION A(LDA,1),WORK(1)
C
C     DCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE
C     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE
C     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK
C     OF A POSITIVE SEMIDEFINITE MATRIX.
C
C     ON ENTRY
C
C         A      DOUBLE PRECISION(LDA,P).
C                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO
C                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED.
C                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED.
C
C         LDA    INTEGER.
C                LDA IS THE LEADING DIMENSION OF THE ARRAY A.
C
C         P      INTEGER.
C                P IS THE ORDER OF THE MATRIX.
C
C         WORK   DOUBLE PRECISION.
C                WORK IS A WORK ARRAY.
C
C         JPVT   INTEGER(P).
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
C                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED.
C                EACH DIAGONAL ELEMENT A(K,K)
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
C                VALUE OF JPVT(K).
C
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
C                                      ELEMENT.
C
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT.
C
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT.
C
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS
C                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO
C                THE BEGINNING OF THE ARRAY A AND FINAL
C                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
C                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE
C                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT
C                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT
C                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.
C                IF JOB .NE. 0, PIVOTING IS DONE.
C
C     ON RETURN
C
C         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR
C                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING.
C
C         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT
C                OF A THAT WAS MOVED INTO THE J-TH POSITION,
C                PROVIDED PIVOTING WAS REQUESTED.
C
C         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL
C                ELEMENT OF THE CHOLESKY FACTOR.
C
C     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN.
C     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL
C     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN
C     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO
C     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE
C     INFO TO BE LESS THAN P.
C
C     LINPACK. THIS VERSION DATED 03/19/79 .
C     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND
C     UNIVERSITY OF MARYLAND.
C
C
C     BLAS DAXPY,DSWAP
C     FORTRAN DSQRT
C
C     INTERNAL VARIABLES
C
      INTEGER PU,PL,PLP1,J,JP,JT,K,KB,KM1,KP1,L,MAXL
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION MAXDIA
      LOGICAL SWAPK,NEGK
C
      PL = 1
      PU = 0
      INFO = P
      IF (JOB .EQ. 0) GO TO 160
C
C        PIVOTING HAS BEEN REQUESTED. REARRANGE THE
C        THE ELEMENTS ACCORDING TO JPVT.
C
         DO 70 K = 1, P
            SWAPK = JPVT(K) .GT. 0
            NEGK = JPVT(K) .LT. 0
            JPVT(K) = K
            IF (NEGK) JPVT(K) = -JPVT(K)
            IF (.NOT.SWAPK) GO TO 60
               IF (K .EQ. PL) GO TO 50
                  CALL DSWAP(PL-1,A(1,K),1,A(1,PL),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PL,PL)
                  A(PL,PL) = TEMP
                  PLP1 = PL + 1
                  IF (P .LT. PLP1) GO TO 40
                  DO 30 J = PLP1, P
                     IF (J .GE. K) GO TO 10
                        TEMP = A(PL,J)
                        A(PL,J) = A(J,K)
                        A(J,K) = TEMP
                     GO TO 20
   10                CONTINUE
                     IF (J .EQ. K) GO TO 20
                        TEMP = A(K,J)
                        A(K,J) = A(PL,J)
                        A(PL,J) = TEMP
   20                CONTINUE
   30             CONTINUE
   40             CONTINUE
                  JPVT(K) = JPVT(PL)
                  JPVT(PL) = K
   50          CONTINUE
               PL = PL + 1
   60       CONTINUE
   70    CONTINUE
         PU = P
         IF (P .LT. PL) GO TO 150
         DO 140 KB = PL, P
            K = P - KB + PL
            IF (JPVT(K) .GE. 0) GO TO 130
               JPVT(K) = -JPVT(K)
               IF (PU .EQ. K) GO TO 120
                  CALL DSWAP(K-1,A(1,K),1,A(1,PU),1)
                  TEMP = A(K,K)
                  A(K,K) = A(PU,PU)
                  A(PU,PU) = TEMP
                  KP1 = K + 1
                  IF (P .LT. KP1) GO TO 110
                  DO 100 J = KP1, P
                     IF (J .GE. PU) GO TO 80
                        TEMP = A(K,J)
                        A(K,J) = A(J,PU)
                        A(J,PU) = TEMP
                     GO TO 90
   80                CONTINUE
                     IF (J .EQ. PU) GO TO 90
                        TEMP = A(K,J)
                        A(K,J) = A(PU,J)
                        A(PU,J) = TEMP
   90                CONTINUE
  100             CONTINUE
  110             CONTINUE
                  JT = JPVT(K)
                  JPVT(K) = JPVT(PU)
                  JPVT(PU) = JT
  120          CONTINUE
               PU = PU - 1
  130       CONTINUE
  140    CONTINUE
  150    CONTINUE
  160 CONTINUE
      DO 270 K = 1, P
C
C        REDUCTION LOOP.
C
         MAXDIA = A(K,K)
         KP1 = K + 1
         MAXL = K
C
C        DETERMINE THE PIVOT ELEMENT.
C
         IF (K .LT. PL .OR. K .GE. PU) GO TO 190
            DO 180 L = KP1, PU
               IF (A(L,L) .LE. MAXDIA) GO TO 170
                  MAXDIA = A(L,L)
                  MAXL = L
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
C
C        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.
C
         IF (MAXDIA .GT. 0.0D0) GO TO 200
            INFO = K - 1
C     ......EXIT
            GO TO 280
  200    CONTINUE
         IF (K .EQ. MAXL) GO TO 210
C
C           START THE PIVOTING AND UPDATE JPVT.
C
            KM1 = K - 1
            CALL DSWAP(KM1,A(1,K),1,A(1,MAXL),1)
            A(MAXL,MAXL) = A(K,K)
            A(K,K) = MAXDIA
            JP = JPVT(MAXL)
            JPVT(MAXL) = JPVT(K)
            JPVT(K) = JP
  210    CONTINUE
C
C        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.
C
         WORK(K) = DSQRT(A(K,K))
         A(K,K) = WORK(K)
         IF (P .LT. KP1) GO TO 260
         DO 250 J = KP1, P
            IF (K .EQ. MAXL) GO TO 240
               IF (J .GE. MAXL) GO TO 220
                  TEMP = A(K,J)
                  A(K,J) = A(J,MAXL)
                  A(J,MAXL) = TEMP
               GO TO 230
  220          CONTINUE
               IF (J .EQ. MAXL) GO TO 230
                  TEMP = A(K,J)
                  A(K,J) = A(MAXL,J)
                  A(MAXL,J) = TEMP
  230          CONTINUE
  240       CONTINUE
            A(K,J) = A(K,J)/WORK(K)
            WORK(J) = A(K,J)
            TEMP = -A(K,J)
            CALL DAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)
  250    CONTINUE
  260    CONTINUE
  270 CONTINUE
  280 CONTINUE
      RETURN
      END

