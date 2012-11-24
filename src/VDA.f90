!
!     Notes to users:  The main program is at the bottom of this file.
!     The code there perform vertex discriminant analysis on the zoo
!     data from the UC Irvine data repository.  VDA.f90 is a self contained
!     and fairly well documented Fortran 90 program.  To run it, you will
!     have to place the file zoo.dat in the same directory as your
!     executable.
!
      MODULE MM_DISCRIMINANT
!
!     Determine double precision and set constants.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBLE = KIND(0.0D0)
      REAL(KIND=DBLE), PARAMETER :: ZERO  = 0.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: ONE   = 1.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: TWO   = 2.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: THREE = 3.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: FOUR  = 4.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: FIVE  = 5.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: SIX   = 6.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: SEVEN = 7.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: EIGHT = 8.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: NINE  = 9.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: TEN   = 10.0_DBLE
      REAL(KIND=DBLE), PARAMETER :: HALF  = ONE/TWO
      REAL(KIND=DBLE), PARAMETER :: PI    = 3.14159265358979_DBLE
!
!     Declare global scalar variables.
!
      REAL(KIND=DBLE) :: DELTA,EPSILON
!
      CONTAINS
!
      SUBROUTINE READ_DATA(CLASS,FEATURE,CASES,FEATURES,INPUT_FILE)
!
!     This subroutine reads in the data and initializes constants and arrays.
!
      IMPLICIT NONE
      INTEGER :: CASES,FEATURES
      INTEGER, ALLOCATABLE, DIMENSION(:) :: CLASS
      REAL(KIND=DBLE), ALLOCATABLE, DIMENSION(:,:) :: COEFFICIENT,FEATURE
      CHARACTER(LEN=100) :: INPUT_FILE
!
      INTEGER :: I,INPUT_UNIT=1
!
!     Read the data.  Fill in the intercept coefficient for each case.
!
      OPEN(UNIT=INPUT_UNIT,FILE=INPUT_FILE)
      DO I = 1,CASES
         READ(INPUT_UNIT,*) FEATURE(I,2:FEATURES+1),CLASS(I)
         FEATURE(I,1) = ONE
      END DO
      CLOSE(INPUT_UNIT)
      END SUBROUTINE READ_DATA

      SUBROUTINE STANDARDIZE_VECTOR(V)
!
!     This subroutine standardizes the vector V so that it has mean 0
!     and variance 1.
!
      IMPLICIT NONE
      REAL(KIND=DBLE) :: CASES,MEAN,STD_DEV
      REAL(KIND=DBLE), DIMENSION(:) :: V
!
      CASES = REAL(SIZE(V),KIND=DBLE)
      MEAN = SUM(V)/CASES
      IF (CASES<=ONE) THEN
         V = ZERO
      ELSE
         STD_DEV = SQRT(SUM((V-MEAN)**2)/(CASES-ONE))
         IF (STD_DEV>ZERO) V = (V-MEAN)/STD_DEV
      END IF
      END SUBROUTINE STANDARDIZE_VECTOR

      SUBROUTINE SYMMETRIZE(MATRIX,UPPER_TO_LOWER)
!
!     This subroutine copies the upper triangle of a matrix onto the
!     lower triangle if the variable UPPER_TO_LOWER is true.  It does
!     the reverse when UPPER_TO_LOWER is false.
!
      IMPLICIT NONE
      INTEGER :: I
      LOGICAL :: UPPER_TO_LOWER
      REAL(KIND=DBLE), DIMENSION(:,:) :: MATRIX
!
      DO I = 2,SIZE(MATRIX,1)
         IF (UPPER_TO_LOWER) THEN
            MATRIX(I,1:I-1) = MATRIX(1:I-1,I)
         ELSE
            MATRIX(1:I-1,I) = MATRIX(I,1:I-1)
         END IF
      END DO
      END SUBROUTINE SYMMETRIZE

      SUBROUTINE FORWARD_SUBSTITUTION(L,A,ERROR)
!
!     This subroutine solves the linear system L X = A by forward
!     substitution following algorithm 3.1.3 of Golub and Van Loan
!     (1996).  The matrix L is lower triangular.  The vector A
!     returns with the solution X.
!
      IMPLICIT NONE
      INTEGER :: J,N
      LOGICAL :: ERROR
      REAL(KIND=DBLE), DIMENSION(:,:) :: L
      REAL(KIND=DBLE), DIMENSION(:) :: A
!
      ERROR = .FALSE.
      N = SIZE(A)
      DO J = 1,N-1
         IF (ABS(L(J,J))<=ZERO) THEN
            ERROR = .TRUE.
            RETURN
         ELSE
            A(J) = A(J)/L(J,J)
            A(J+1:N) = A(J+1:N)-A(J)*L(J+1:N,J)
         END IF
      END DO
      IF (ABS(L(N,N))<=ZERO) THEN
         ERROR = .TRUE.
      ELSE
         A(N) = A(N)/L(N,N)
      END IF
      END SUBROUTINE FORWARD_SUBSTITUTION

      SUBROUTINE BACKWARD_SUBSTITUTION(U,B,ERROR)
!
!     This subroutine solves the linear system U X = B by backward
!     substitution following algorithm 3.1.4 of Golub and Van Loan
!     (1996).  The matrix U is upper triangular.  The vector B
!     returns with the solution X.
!
      IMPLICIT NONE
      INTEGER :: J,N
      LOGICAL :: ERROR
      REAL(KIND=DBLE), DIMENSION(:,:) :: U
      REAL(KIND=DBLE), DIMENSION(:) :: B
!
      ERROR = .FALSE.
      N = SIZE(U,2)
      DO J = N,2,-1
         IF (ABS(U(J,J))<=ZERO) THEN
            ERROR = .TRUE.
            RETURN
         ELSE
            B(J) = B(J)/U(J,J)
            B(1:J-1) = B(1:J-1)-B(J)*U(1:J-1,J)
         END IF
      END DO
      IF (ABS(U(1,1))<=ZERO) THEN
         ERROR = .TRUE.
      ELSE
         B(1) = B(1)/U(1,1)
      END IF
      END SUBROUTINE BACKWARD_SUBSTITUTION

      SUBROUTINE CHOLESKY(MATRIX,ERROR)
!
!     This subroutine overwrites the lower triangle of a positive
!     definite matrix with its Cholesky decomposition following
!     algorithm 4.2.1 of Golub and Van Loan (1966).  If a nonpositive
!     diagonal entry is encountered at any stage, then the subroutine
!     aborts with an error.
!
      IMPLICIT NONE
      INTEGER :: J,K,N
      LOGICAL :: ERROR
      REAL(KIND=DBLE) :: T
      REAL(KIND=DBLE), DIMENSION(:,:) :: MATRIX
!
      ERROR = .FALSE.
      N = SIZE(MATRIX,1)
      DO J = 1,N
         DO K = 1,J-1
            MATRIX(J:N,J) = MATRIX(J:N,J)-MATRIX(J:N,K)*MATRIX(J,K)
         END DO
         IF (MATRIX(J,J)>ZERO) THEN
            MATRIX(J:N,J) = MATRIX(J:N,J)/SQRT(MATRIX(J,J))
         ELSE
            ERROR = .TRUE.
            RETURN
         END IF
      END DO
      END SUBROUTINE CHOLESKY

      SUBROUTINE VERTEX_DISCRIMINANT_ANALYSIS(CLASS,COEFFICIENT,FEATURE, &
         PREDICTED,LAMBDA,TRAINING_ERROR_RATE,CASES,CLASSES,FEATURES,ERROR)
!
!     This subroutine carries out vertex discriminant analysis by
!     the MM algorithm.  The user is responsible for setting the
!     tuning constant LAMBDA, the number of CASES, the number of
!     CLASSES, the number of FEATURES, the FEATURE matrix, and
!     the CLASS vector.  On exit the subroutine returns the
!     TRAINING_ERROR_RATE, the PREDICTED class of each case, and
!     the estimated COEFFICIENT matrix for classifying new cases.
!     The rows of FEATURE correspond to cases and the columns to
!     features.  All entries of the first column should equal 1.
!     These entries serve in estimating the intercept parameters.
!     The entries of CLASS should be integers between 1 and CLASSES.
!     Users wanting to classify new cases should reconstitute the
!     vertices of the regular simplex using the code below.  The
!     first column of COEFFICIENT gives the estimated intercepts.
!     For more details see the paper:  Lange K, Wu TT (2007) An MM
!     algorithm for multicategory vertex discriminant analysis. J
!     Computational Graphical Stat.
!
      IMPLICIT NONE
!
!     Declare scalars in the calling sequence.
!
      INTEGER :: CASES,CLASSES,FEATURES
      LOGICAL :: ERROR
      REAL(KIND=DBLE) :: LAMBDA,TRAINING_ERROR_RATE
!
!     Declare arrays in the calling sequence.
!
      INTEGER, DIMENSION(CASES) :: CLASS,PREDICTED
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,FEATURES+1) :: COEFFICIENT
      REAL(KIND=DBLE), DIMENSION(CASES,FEATURES+1) :: FEATURE
!
!     Declare local scalars.
!
      INTEGER :: BEST_CLASS,I,J,K,ITERATION,MAX_ITERATIONS = 5000
      REAL(KIND=DBLE) :: A,B,C,D,L2_NORM,LOSS,OLD_LOSS,PENALTY,W
      REAL(KIND=DBLE) :: BEST_NORM_SQUARED,NORM_SQUARED,PENALTY_CONSTANT
      REAL(KIND=DBLE) :: CONVERGENCE_CRITERION = TEN**(-5),EPSILON
      REAL(KIND=DBLE) :: DELTA = TEN**(-5),EXPAND = TWO
!
!     Declare local arrays.
!
      REAL(KIND=DBLE), DIMENSION(CLASSES,CLASSES-1) :: VERTEX
      REAL(KIND=DBLE), DIMENSION(CLASSES-1) :: CENTER,PREDICTOR,RESIDUAL
      REAL(KIND=DBLE), DIMENSION(FEATURES+1) :: ESTIMATE
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,FEATURES+1) :: XTWY
      REAL(KIND=DBLE), DIMENSION(FEATURES+1,FEATURES+1) :: XTWX
!
!     Construct the vertices of a regular simplex in a Euclidean space with
!     dimension equal to CLASSES-1.
!
      A = REAL(CLASSES,KIND=DBLE)
      B = A-ONE
      C = SQRT(A)
      D = SQRT(B)
      VERTEX(1,:) = ONE/D
      DO I = 1,CLASSES-1
         VERTEX(I+1,:) = -(ONE+C)/(B*D)
         VERTEX(I+1,I) = VERTEX(I+1,I)+C/D
      END DO
      EPSILON = HALF*SQRT(TWO*REAL(CLASSES,KIND=DBLE)/(REAL(CLASSES,KIND=DBLE)-ONE))
!
!     Initialize the parameters, the previous value of the loss function, and 
!     the penalty constant.
!
      COEFFICIENT = ZERO
      OLD_LOSS = HUGE(ONE)
      PENALTY_CONSTANT = LAMBDA*REAL(CASES,KIND=DBLE)
!
!     Enter the MM iteration loop.
!
      DO ITERATION = 1,MAX_ITERATIONS
!
!     Initialize variables.
!
         LOSS = ZERO
         XTWX = ZERO
         XTWY = ZERO
!
!     Loop over all cases.
!
         DO I = 1,CASES
!
!     Compute the current residual and its contribution to the epsilon insensitive loss.
!
            RESIDUAL = VERTEX(CLASS(I),:)-MATMUL(COEFFICIENT,FEATURE(I,:))
            L2_NORM = SQRT(SUM(RESIDUAL**2))
            LOSS = LOSS+MAX(L2_NORM-EPSILON,ZERO)
!
!     Compute the center and weight of the quadratic surrogate for the loss.
!
            IF (L2_NORM>=TWO*EPSILON) THEN
               CENTER = ZERO
               W = TWO*L2_NORM
            ELSE IF (L2_NORM<=EPSILON) THEN
               CENTER = RESIDUAL
               W = FOUR*MAX(EPSILON-L2_NORM,DELTA)
            ELSE
               CENTER = (TWO*EPSILON/L2_NORM-ONE)*RESIDUAL
               W = FOUR*MAX(L2_NORM-EPSILON,DELTA)
            END IF
            W = ONE/W
!
!     Update XTWX and XTWY.  These are the standard entities in weighted
!     least squares.
!
            DO J = 1,FEATURES+1
               A = FEATURE(I,J)
               DO K = 1,J
                  XTWX(J,K) = XTWX(J,K)+A*W*FEATURE(I,K)
               END DO
               XTWY(:,J) = XTWY(:,J)+A*W*(VERTEX(CLASS(I),:)-CENTER)
            END DO
         END DO
!
!     Add the penalty terms to XTWX.
!
         DO I = 2,FEATURES+1
            XTWX(I,I) = XTWX(I,I)+PENALTY_CONSTANT
         END DO
!
!     Add the penalty terms to the loss function.
!
         PENALTY = ZERO
         DO I = 1,CLASSES-1
            A = DOT_PRODUCT(COEFFICIENT(I,2:FEATURES+1),COEFFICIENT(I,2:FEATURES+1))
            PENALTY = PENALTY+A
         END DO
         LOSS = LOSS+PENALTY_CONSTANT*PENALTY
!
!     Find the Cholesky decomposition of XTWX.
!
         CALL CHOLESKY(XTWX,ERROR)
         IF (ERROR) EXIT
         CALL SYMMETRIZE(XTWX,UPPER_TO_LOWER=.FALSE.)
!
!     Solve the separated least squares problems.
!
         DO I = 1,CLASSES-1
!
!     Use forward and backward substitution.
!
            ESTIMATE = XTWY(I,:)
            CALL FORWARD_SUBSTITUTION(XTWX,ESTIMATE,ERROR)
            IF (ERROR) EXIT
            CALL BACKWARD_SUBSTITUTION(XTWX,ESTIMATE,ERROR)
            IF (ERROR) EXIT
!
!     EXPAND = 2 corresponds to step doubling.
!
            IF (ITERATION<=10) THEN
               COEFFICIENT(I,:) = ESTIMATE
            ELSE
               COEFFICIENT(I,:) = COEFFICIENT(I,:)+EXPAND*(ESTIMATE-COEFFICIENT(I,:))
            END IF
         END DO
!
!     Check for convergence.
!
         IF (ABS(OLD_LOSS-LOSS)<CONVERGENCE_CRITERION) THEN
            EXIT
         ELSE
            OLD_LOSS = LOSS
         END IF
!
      END DO
!
!     Check for a search error.
!
      IF (ERROR) THEN
         TRAINING_ERROR_RATE = HUGE(ONE)
!
!     Reclassify cases and compute the TRAINING_ERROR_RATE.
!
      ELSE
         TRAINING_ERROR_RATE = ZERO
         DO I = 1,CASES
            BEST_NORM_SQUARED = HUGE(ONE)
            PREDICTOR = MATMUL(COEFFICIENT,FEATURE(I,:))
            DO J = 1,CLASSES
               NORM_SQUARED = SUM((VERTEX(J,:)-PREDICTOR)**2)
               IF (NORM_SQUARED<BEST_NORM_SQUARED) THEN
                  BEST_CLASS = J
                  BEST_NORM_SQUARED = NORM_SQUARED
               END IF
            END DO
            IF (CLASS(I)/=BEST_CLASS) THEN
               TRAINING_ERROR_RATE = TRAINING_ERROR_RATE+ONE
            END IF
            PREDICTED(I) = BEST_CLASS
         END DO
         TRAINING_ERROR_RATE = TRAINING_ERROR_RATE/REAL(CASES,KIND=DBLE)
      END IF
      END SUBROUTINE VERTEX_DISCRIMINANT_ANALYSIS

      SUBROUTINE FETCH_OBJECTIVE(FEATURE,CASES,CLASSES,FEATURES,LAMBDA1,LAMBDA2,COEFFICIENT,&
           RESIDUAL,D1OBJECTIVE,D2OBJECTIVE,OBJECTIVE,I,J,FIRST)
!
!     Fetch the objective function and its first two derivatives.  I is the
!     dimension in vertex space and J is the feature.
!
      IMPLICIT NONE
!
!     Declare scalars in the calling sequence.
!
      INTEGER :: I,J,CASES,CLASSES,FEATURES
      LOGICAL, OPTIONAL :: FIRST
      REAL(KIND=DBLE) :: D1OBJECTIVE,D2OBJECTIVE,OBJECTIVE,LAMBDA1,LAMBDA2
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,FEATURES+1) :: COEFFICIENT
      REAL(KIND=DBLE), DIMENSION(CASES,FEATURES+1) :: FEATURE
!
!     Declare local scalar variables.
!
      INTEGER :: K
      REAL(KIND=DBLE), SAVE :: C1,C2,C3,C4,C5,MU1,MU2
      REAL(KIND=DBLE) :: D1P,D2P,NORM,NORMSQ,P,RX,RX2,S,S2,S3,T,X
!
!     Declare local arrays.
!
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,CASES) :: RESIDUAL

!
!     On the first call, set constants pertinent to all remaining calls.
!
      IF (PRESENT(FIRST)) THEN
         C1 = EPSILON-DELTA
         C2 = EPSILON+DELTA
         C3 = ONE/(FOUR**2*DELTA**3)
         C4 = ONE/(FOUR*DELTA**3)
         C5 = THREE/(FOUR*DELTA**3)
         MU1 = CASES*LAMBDA1
         MU2 = CASES*LAMBDA2
         !RETURN
      END IF
!
!     Initialize the objective function.
!
      OBJECTIVE = ZERO
      D1OBJECTIVE = ZERO
      D2OBJECTIVE = ZERO
!
!     Loop over all cases K.
!
      DO K = 1,CASES
         NORMSQ = SUM(RESIDUAL(:,K)**2)
         NORM = SQRT(NORMSQ)
!
!     Skip the loss for the current case if its predicted value is within
!     a distance EPSILON-DELTA of its associated vertex.
!
         IF (NORM<C1) CYCLE
         X = FEATURE(K,J)
         RX = RESIDUAL(I,K)*X
!
!     Compute the loss for the current case if its predicted value is beyond
!     a distance EPSILON+DELTA of its associated vertex.
!
         IF (NORM>C2) THEN
            OBJECTIVE = OBJECTIVE+NORM-EPSILON
            S = ONE/NORM
            D1OBJECTIVE = D1OBJECTIVE-S*RX
            D2OBJECTIVE = D2OBJECTIVE+S*(X*X-RX*RX/NORMSQ)
!
!     Compute the loss for the current case if its predicted distance falls in
!     the transition zone for its associated vertex.
!
         ELSE
            S = NORM-C1
            T = C2-NORM
            D2P = C5*S*T
            S2 = S*S
            T = T+DELTA
            D1P = C4*S2*T
            T = T+DELTA
            S3 = S2*S
            P = C3*S3*T
            OBJECTIVE = OBJECTIVE+P
            S = D1P/NORM
            D1OBJECTIVE = D1OBJECTIVE-S*RX
            T = D2P/NORMSQ
            RX2 = RX*RX
            D2OBJECTIVE = D2OBJECTIVE+T*RX2+S*(X*X-RX2/NORMSQ)
         END IF
      END DO
!
!     Compute the PENALTY2 for the current parameter.
!
      IF (J>1) THEN
         NORMSQ = SUM(COEFFICIENT(:,J)**2)
         NORM = SQRT(NORMSQ)
         OBJECTIVE = OBJECTIVE+MU2*NORM+MU1*SUM(ABS(COEFFICIENT(:,J)))
!
!     If the parameter is parked at 0, then reset the first derivative of the objective
!     so that it either remains at 0 or moves in the correct direction.
!
         IF (NORM<=ZERO) THEN
            IF (ABS(D1OBJECTIVE)<=MU2+MU1) THEN
               D1OBJECTIVE = ZERO
            ELSE IF (D1OBJECTIVE+MU2+MU1<=ZERO) THEN
               D1OBJECTIVE = D1OBJECTIVE+MU2+MU1
            ELSE
               D1OBJECTIVE = D1OBJECTIVE-MU2-MU1
            END IF
!
!     Deal with the Euclidean PENALTY2 away from the origin.
!
         ELSE
            S = COEFFICIENT(I,J)
            T = S/NORM
            D1OBJECTIVE = D1OBJECTIVE+MU2*T
            D2OBJECTIVE = D2OBJECTIVE+MU2*(ONE-T*T)/NORM
            IF (ABS(D1OBJECTIVE)<=MU1) THEN
               D1OBJECTIVE = ZERO
            ELSE IF (D1OBJECTIVE+MU1<=ZERO) THEN
               D1OBJECTIVE = D1OBJECTIVE+MU1
            ELSE
               D1OBJECTIVE = D1OBJECTIVE-MU1
            END IF
         END IF
      END IF
      END SUBROUTINE FETCH_OBJECTIVE

     SUBROUTINE CYCLIC_EUCLIDEAN_VDA(FEATURE,CLASS,CASES,CLASSES,FEATURES,LAMBDA1,&
          LAMBDA2,POST_CLASS,COEFFICIENT,TRAINING_ERROR_RATE,ERROR)
!
!     This subroutine carries out cyclic coordinate descent in the context
!     of vertex discriminant analysis with a Euclidean PENALTY2.  The user
!     is responsible for setting the PENALTY2 constant LAMBDA2.  At the top
!     of the data file, define the number of CASES, the number of CLASSES,
!     and the number of FEATURES.  Then list the FEATURE matrix, and the
!     CLASS vector.  On exit the current subroutine returns the predicted
!     class POST_CLASS of each case, the estimated COEFFICIENT matrix for
!     classifying new cases, and the TRAINING_ERROR_RATE.  The rows of
!     FEATURE correspond to cases and the columns to features.  All entries
!     of the first column should equal 1.  These entries serve in estimating
!     the intercept parameters.  The entries of CLASS should be integers
!     between 1 and CLASSES.  Users wanting to classify new cases should
!     reconstitute the vertices of the regular simplex using the code below.
!
      IMPLICIT NONE
!
!     Declare scalars in the calling sequence.
!
      LOGICAL :: ERROR
      REAL(KIND=DBLE) :: TRAINING_ERROR_RATE,LAMBDA1,LAMBDA2
      INTEGER :: CASES,CLASSES,FEATURES
      INTEGER, DIMENSION(CASES) :: CLASS,POST_CLASS
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,FEATURES+1) :: COEFFICIENT
      REAL(KIND=DBLE), DIMENSION(CASES,FEATURES+1) :: FEATURE
!
!     Declare local scalars.
!
      INTEGER :: BEST_CLASS,I,ITERATION,J,K,MAX_ITERATIONS
      INTEGER :: MAX_STEPS,NEWTON_ITERATION,STEP
      REAL(KIND=DBLE) :: FULL_OBJECTIVE,CRITERION,LOSS,MU1,MU2,PENALTY2,PENALTY1
      REAL(KIND=DBLE) :: D1OBJECTIVE,D2OBJECTIVE,NEW_OBJECTIVE,OBJECTIVE
      REAL(KIND=DBLE) :: BEST_NORMSQ,NORM,NORMSQ,TEST
      REAL(KIND=DBLE) :: A,B,C,D,D1P,D2P,P
!
!     Declare local arrays.
!
      
      REAL(KIND=DBLE), DIMENSION(CLASSES-1) :: PREDICTED
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,CASES) :: RESIDUAL
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,CLASSES) :: VERTEX
!
!     Construct the vertices of a regular simplex in a Euclidean space with
!     dimension equal to CLASSES-1.
!
      A = REAL(CLASSES,KIND=DBLE)
      B = A-ONE
      C = SQRT(A)
      D = SQRT(B)
      VERTEX(:,1) = ONE/D
      DO I = 1,CLASSES-1
         VERTEX(:,I+1) = -(ONE+C)/(B*D)
         VERTEX(I,I+1) = VERTEX(I,I+1)+C/D
      END DO
!
!     Define EPSILON and DELTA used in the definition of smoothed epsilon
!     insensitive distance. Also define the convergence CRITERION, the
!     adjusted PENALTY2 constant MU2, the maximum number of outer iterations,
!     and the maximum number of step halvings for the backtrack loop.  Also
!     initialize the error flag and the full objective function.
!
      EPSILON = HALF*SQRT(TWO*A/(A-ONE))
      DELTA = ONE/(TWO*TEN)
      CRITERION = TEN**(-4)
      MU1 = CASES*LAMBDA1
      MU2 = CASES*LAMBDA2
      MAX_ITERATIONS = 1000
      MAX_STEPS = 5
      ERROR = .FALSE.
      FULL_OBJECTIVE = HUGE(ONE)
!
!     Initialize the parameters and the residuals.
!
      COEFFICIENT = ZERO
      DO I = 1,CASES
         RESIDUAL(:,I) = VERTEX(:,CLASS(I))
      END DO
!
!     Initialize the constants in the objective function routine.
!
      CALL FETCH_OBJECTIVE(FEATURE,CASES,CLASSES,FEATURES,LAMBDA1,&
           LAMBDA2,COEFFICIENT,RESIDUAL,A,B,C,1,1,.TRUE.)
      PENALTY1 = ZERO
      PENALTY2 = ZERO
      DO J = 2,FEATURES+1
         PENALTY1 = PENALTY1+MU1*SUM(ABS(COEFFICIENT(:,J)))
         PENALTY2 = PENALTY2+MU2*SQRT(SUM(COEFFICIENT(:,J)**2))
      END DO
      !WRITE(*,"(I5,1X,5(F10.3,1X)/)") 0,C+PENALTY2+PENALTY1,C,PENALTY1,PENALTY2
!
!     Enter the outer iteration loop.
!
      DO ITERATION = 1,MAX_ITERATIONS
!
!     Update the parameters cyclically.
!
         FEATURE_LOOP: DO J = 1,FEATURES+1
            CLASS_LOOP: DO I = 1,CLASSES-1
!
!     Update the current parameter by Newton's method.  Skip a parameter
!     parked at 0 if there is insufficient push.
!
               CALL FETCH_OBJECTIVE(FEATURE,CASES,CLASSES,FEATURES,LAMBDA1, LAMBDA2,&
                    COEFFICIENT,RESIDUAL,D1OBJECTIVE,D2OBJECTIVE,OBJECTIVE,I,J)
               IF (ABS(D1OBJECTIVE)<=ZERO) CYCLE CLASS_LOOP
!
!     Enter the Newton loop.
!
               NEWTON_LOOP: DO NEWTON_ITERATION = 1,10
!
!     P is the current value of the parameter, and D is the Newton increment.
!
                  P = COEFFICIENT(I,J)
                  D = -D1OBJECTIVE/D2OBJECTIVE
!
!     Compute a new function value.  If it represents a sufficient decrease,
!     then accept the new point.  Otherwise step halve.
!
                  BACKTRACK_LOOP: DO STEP = 0,MAX_STEPS
                     A = P+D-COEFFICIENT(I,J)
                     COEFFICIENT(I,J) = P+D
!
!     Update the residual vector.
!
                     DO K = 1,CASES
                        RESIDUAL(I,K) = RESIDUAL(I,K)-A*FEATURE(K,J)
                     END DO
!
!     Compute a new value of the objective function and check whether it
!     gives an improvement.  If not, step halve.
!
                     CALL FETCH_OBJECTIVE(FEATURE,CASES,CLASSES,FEATURES,LAMBDA1,&
                          LAMBDA2,COEFFICIENT,RESIDUAL,D1OBJECTIVE,D2OBJECTIVE, &
                        NEW_OBJECTIVE,I,J)
!
!    objective is wrong here; PENALTY2 is miscomputed.
!
                     IF (NEW_OBJECTIVE<=OBJECTIVE) EXIT BACKTRACK_LOOP
                     D = HALF*D
                  END DO BACKTRACK_LOOP
!
!     Test for convergence.
!
                  TEST = ABS(OBJECTIVE-NEW_OBJECTIVE)
                  OBJECTIVE = NEW_OBJECTIVE
                  IF (TEST<=TEN**(-8)) EXIT
               END DO NEWTON_LOOP
!
!     Reset parameters close to 0 to 0.
!
               IF (ABS(COEFFICIENT(I,J))<=TEN**(-8)) COEFFICIENT(I,J) = ZERO
            END DO CLASS_LOOP
         END DO FEATURE_LOOP
!
!     Compute the loss function.  Note that the objective function delivers
!     the PENALTY2 for only the current feature.  For the first feature there
!     is no PENALTY2.  Hence, compute the full PENALTY2.
!
         CALL FETCH_OBJECTIVE(FEATURE,CASES,CLASSES,FEATURES,LAMBDA1,&
              LAMBDA2,COEFFICIENT,RESIDUAL,D1OBJECTIVE,D2OBJECTIVE,LOSS,1,1)
         PENALTY1 = ZERO
         PENALTY2 = ZERO
         DO J = 2,FEATURES+1
            PENALTY1 = PENALTY1+MU1*SUM(ABS(COEFFICIENT(:,J)))
            PENALTY2 = PENALTY2+MU2*SQRT(SUM(COEFFICIENT(:,J)**2))
         END DO
!
!     Output the iteration number and the value of the objective function.
!
         OBJECTIVE = LOSS+PENALTY2+PENALTY1
         !IF (ITERATION==1.OR.MOD(ITERATION,10)==0) THEN
            !PRINT*," ITERATION = ",ITERATION," FUN = ",OBJECTIVE
            !WRITE(*,"(I5,1X,5(F10.3,1X)/)") ITERATION,OBJECTIVE,LOSS,PENALTY1,PENALTY2
         !END IF
!
!     Check for a descent failure.
!
         IF (OBJECTIVE>FULL_OBJECTIVE) THEN
            ERROR = .TRUE.
            RETURN
            PRINT*," *** ERROR *** OBJECTIVE FUNCTION INCREASE"
            EXIT
         END IF
!
!     Check for convergence.  If convergence is not declared, then reset
!     the full objective function.
!
         IF (FULL_OBJECTIVE-OBJECTIVE<CRITERION) EXIT
         FULL_OBJECTIVE = OBJECTIVE
      END DO
!
!     Reclassify cases and compute the TRAINING_ERROR_RATE.
!
      TRAINING_ERROR_RATE = ZERO
      DO I = 1,CASES
         BEST_NORMSQ = HUGE(ONE)
         PREDICTED = VERTEX(:,CLASS(I))-RESIDUAL(:,I)
         DO J = 1,CLASSES
            NORMSQ = SUM((VERTEX(:,J)-PREDICTED)**2)
            IF (NORMSQ<BEST_NORMSQ) THEN
               BEST_CLASS = J
               BEST_NORMSQ = NORMSQ
            END IF
         END DO
         IF (CLASS(I)/=BEST_CLASS) THEN
            TRAINING_ERROR_RATE = TRAINING_ERROR_RATE+ONE
         END IF
         POST_CLASS(I) = BEST_CLASS
      END DO
      TRAINING_ERROR_RATE = TRAINING_ERROR_RATE/REAL(CASES,KIND=DBLE)

      END SUBROUTINE CYCLIC_EUCLIDEAN_VDA

      END MODULE MM_DISCRIMINANT

   
      SUBROUTINE  VDA(FEATURE,CLASS,CASES,CLASSES,FEATURES, &
           LAMBDA,PREDICTED,COEFFICIENT,TRAINING_ERROR_RATE)

      USE MM_DISCRIMINANT
      IMPLICIT NONE
      INTEGER :: CASES,CLASSES,FEATURES,I
      LOGICAL :: ERROR
      REAL(KIND=DBLE) :: LAMBDA,TRAINING_ERROR_RATE
!
!     Declare and allocate model arrays.
!
      INTEGER, DIMENSION(CASES) :: CLASS,PREDICTED
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,FEATURES+1) :: COEFFICIENT
      REAL(KIND=DBLE), DIMENSION(CASES,FEATURES+1) :: FEATURE 
!
!     Standarize it.
!
      DO I = 2,FEATURES+1
        CALL STANDARDIZE_VECTOR(FEATURE(:,I))
      END DO
!
!     Set the tuning constant LAMBDA and perform discriminant analysis.
!
      CALL VERTEX_DISCRIMINANT_ANALYSIS(CLASS,COEFFICIENT,FEATURE, &
         PREDICTED,LAMBDA,TRAINING_ERROR_RATE,CASES,CLASSES,FEATURES,ERROR)
!
      END SUBROUTINE VDA

      SUBROUTINE VDA_LE(FEATURE,CLASS,CASES,CLASSES,FEATURES,LAMBDA1,&
           LAMBDA2,POST_CLASS,COEFFICIENT,TRAINING_ERROR_RATE)
!
!     This program performs Vertex Discriminant Analysis on the zoo data
!     from the UC Irvine data repository.  For more documentation, see
!     the preamble to the subroutine CYCLIC_EUCLIDEAN_VDA.
!
      USE MM_DISCRIMINANT
      IMPLICIT NONE
!
!     Declare local scalar variables.
!
      INTEGER :: I,J,NONZEROS,CASES,CLASSES,FEATURES
      LOGICAL :: ERROR
      REAL(KIND=DBLE) :: TRAINING_ERROR_RATE,LAMBDA1,LAMBDA2
      INTEGER, DIMENSION(CASES) :: CLASS,POST_CLASS
      REAL(KIND=DBLE), DIMENSION(CLASSES-1,FEATURES+1) :: COEFFICIENT
      REAL(KIND=DBLE), DIMENSION(CASES,FEATURES+1) :: FEATURE
!
!     Read the data and standarize all features.
!
      DO I = 2,FEATURES+1
        CALL STANDARDIZE_VECTOR(FEATURE(:,I))
      END DO

      CALL CYCLIC_EUCLIDEAN_VDA(FEATURE,CLASS,CASES,CLASSES,FEATURES,LAMBDA1,&
              LAMBDA2,POST_CLASS,COEFFICIENT,TRAINING_ERROR_RATE,ERROR)

      END SUBROUTINE VDA_LE
