      SUBROUTINE DTRSL(T,LDT,N,B,JOB,INFO)
      INTEGER LDT,N,JOB,INFO
      DOUBLE PRECISION T(LDT,1),B(1)
C
C
C     DTRSL SOLVES SYSTEMS OF THE FORM
C
C                   T * X = B
C     OR
C                   TRANS(T) * X = B
C
C     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE TRANS(T)
C     DENOTES THE TRANSPOSE OF THE MATRIX T.
C
C     ON ENTRY
C
C         T         DOUBLE PRECISION(LDT,N)
C                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO
C                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                   USED TO STORE OTHER INFORMATION.
C
C         LDT       INTEGER
C                   LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C
C         N         INTEGER
C                   N IS THE ORDER OF THE SYSTEM.
C
C         B         DOUBLE PRECISION(N).
C                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.
C
C         JOB       INTEGER
C                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.
C                   IF JOB IS
C
C                        00   SOLVE T*X=B, T LOWER TRIANGULAR,
C                        01   SOLVE T*X=B, T UPPER TRIANGULAR,
C                        10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR,
C                        11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR.
C
C     ON RETURN
C
C         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.
C                   OTHERWISE B IS UNALTERED.
C
C         INFO      INTEGER
C                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.
C                   OTHERWISE INFO CONTAINS THE INDEX OF
C                   THE FIRST ZERO DIAGONAL ELEMENT OF T.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C     FORTRAN MOD
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,TEMP
      INTEGER CASE,J,JJ
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 150
C
C        CHECK FOR ZERO DIAGONAL ELEMENTS.
C
         DO 10 INFO = 1, N
C     ......EXIT
            IF (T(INFO,INFO) .EQ. 0.0D0) GO TO 150
   10    CONTINUE
         INFO = 0
C
C        DETERMINE THE TASK AND GO TO IT.
C
         CASE = 1
         IF (MOD(JOB,10) .NE. 0) CASE = 2
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2
         GO TO (20,50,80,110), CASE
C
C        SOLVE T*X=B FOR T LOWER TRIANGULAR
C
   20    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 40
            DO 30 J = 2, N
               TEMP = -B(J-1)
               CALL DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
               B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
         GO TO 140
C
C        SOLVE T*X=B FOR T UPPER TRIANGULAR.
C
   50    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 70
            DO 60 JJ = 2, N
               J = N - JJ + 1
               TEMP = -B(J+1)
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
         GO TO 140
C
C        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.
C
   80    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 100
            DO 90 JJ = 2, N
               J = N - JJ + 1
               B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)
               B(J) = B(J)/T(J,J)
   90       CONTINUE
  100       CONTINUE
         GO TO 140
C
C        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.
C
  110    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 130
            DO 120 J = 2, N
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
               B(J) = B(J)/T(J,J)
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
