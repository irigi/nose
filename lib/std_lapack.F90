!
! Simple interface to lapack based on original LAPACK95 package
! See: LAPACK95 Users' Guide, available at http://www.netlib.org/lapack95/lug95/
!
! The routines std_lapack_geev and std_lapack_gesv are more or less cut and pasted versions of
! la_geev and la_gesv subrotines of that package. The LAPACK95 licence seems to allow this.
!
! std_types are replacing the LAPACK95 types
! subroutine LSAME and ERINFO are used in their original version
!
! Author: Tomas Mancal, mancal@karlov.mff.cuni.cz
! Date: 2008/02/15
!   
module std_lapack
 
    use std_types
    
    IMPLICIT NONE

    ! Interface to lapack routines
    
    INTERFACE LA_GEEV_m

       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )
         use std_types
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(DP), INTENT(INOUT) :: A(LDA,*)
         REAL(DP), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), WR(*), WI(*), &
                                 WORK(*)
      END SUBROUTINE DGEEV

       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,&
                        WORK, LWORK, RWORK, INFO )
                        use std_types
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(DPC), INTENT(OUT) :: RWORK(*)
         COMPLEX(DPC), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(DPC), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), W(*),      &
                                    WORK(*)
      END SUBROUTINE ZGEEV

    END INTERFACE

    INTERFACE LA_GESV_m

       SUBROUTINE DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
                        use std_types
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         REAL(DP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE DGESV

       SUBROUTINE ZGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
                        use std_types
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         COMPLEX(DPC), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE ZGESV

    END INTERFACE

    INTERFACE std_lapack_gesv
    	MODULE PROCEDURE std_lapack_gesv_cmplx
    	MODULE PROCEDURE std_lapack_gesv_real
    END INTERFACE

    INTERFACE std_lapack_geev
    	MODULE PROCEDURE std_lapack_geev_cmplx
    	MODULE PROCEDURE std_lapack_geev_real
    END INTERFACE

contains


!    subroutine std_lapack_geev(AC,W,W2,S1,info)
!        real(dp), dimension(:,:), intent(out) :: S1
!        real(dp), dimension(:) :: W,W2
!        real(dp), dimension(:,:):: AC
!        integer :: info

!        call la_geev_i(AC,W,W2,VR=S1,INFO=info)
!    end subroutine std_lapack_geev    
    
!    subroutine std_lapack_gesv(Ai,B,info)
!        real(dp), dimension(:,:), intent(out) :: Ai
!        real(dp), dimension(:,:) :: B
!        integer :: info
!        call la_gesv_i(Ai,B,INFO=info)
!    end subroutine std_lapack_gesv    


    subroutine std_lapack_geev_cmplx( A, W, W2, VL, VR, INFO )
       IMPLICIT NONE
    !  .. SCALAR ARGUMENTS ..
       INTEGER, INTENT(OUT), OPTIONAL :: INFO
    !  .. ARRAY ARGUMENTS ..
       COMPLEX(dpc), INTENT(INOUT) :: A(:,:)
       COMPLEX(dpc), INTENT(OUT) :: W(:), W2(:)
       COMPLEX(dpc), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)

    !  .. LOCAL PARAMETERS ..
       CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEEV'
    !  .. LOCAL SCALARS ..
       CHARACTER(LEN=1) :: LJOBVL, LJOBVR
       INTEGER, SAVE :: LWORK = 0
       INTEGER :: N, NN, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR
    !  .. LOCAL ARRAYS ..
       COMPLEX(dpc), TARGET :: LLVL(1,1), LLVR(1,1)
       COMPLEX(dpc), POINTER :: WORK(:)
       REAL(dp), POINTER :: RWORK(:) !!!!!!!!!!!!!!!!
    !  .. INTRINSIC FUNCTIONS ..
       INTRINSIC MAX, PRESENT, SIZE
    !  .. EXECUTABLE STATEMENTS ..
       LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
       IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
       ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
       IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
       ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
    !  .. TEST THE ARGUMENTS
       IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
       !HACK
       ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
       !ELSE IF( SIZE( WI ) /= N )THEN; LINFO = -3
       ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -4
       ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -5
       ELSE IF( N > 0 )THEN
          ALLOCATE(RWORK(2*N))
          NN = 3; IF( LSAME(LJOBVL,'V').OR.LSAME(LJOBVR,'V') ) NN = NN + 1
          LWORK = MAX( 1, NN*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
             IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
                LWORK = MAX( 1, NN*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
                IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
             END IF
          IF( ISTAT == 0 ) THEN
            IF( PRESENT(VL) )THEN
               IF( PRESENT(VR) )THEN
                   CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            VL, S1VL, VR, S1VR, WORK, LWORK, RWORK, LINFO )
           ELSE
              CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            VL, S1VL, LLVR, S1VR, WORK, LWORK, RWORK, LINFO )
           ENDIF
         ELSE
           IF( PRESENT(VR) )THEN
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            LLVL, S1VL, VR, S1VR, WORK, LWORK, RWORK, LINFO )
           ELSE
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, W, &
                            LLVL, S1VL, LLVR, S1VR, WORK, LWORK, RWORK, LINFO )
           ENDIF
         ENDIF
             IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
          ELSE; LINFO = -100; ENDIF
          DEALLOCATE(WORK, STAT=ISTAT1)
       ENDIF
       CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
    end subroutine std_lapack_geev_cmplx


    subroutine std_lapack_geev_real( A, WR, WI, VL, VR, INFO )
       IMPLICIT NONE
    !  .. SCALAR ARGUMENTS ..
       INTEGER, INTENT(OUT), OPTIONAL :: INFO
    !  .. ARRAY ARGUMENTS ..
       REAL(dp), INTENT(INOUT) :: A(:,:)
       REAL(dp), INTENT(OUT) :: WR(:), WI(:)
       REAL(dp), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
    
    !  .. LOCAL PARAMETERS ..
       CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEEV'
    !  .. LOCAL SCALARS ..
       CHARACTER(LEN=1) :: LJOBVL, LJOBVR
       INTEGER, SAVE :: LWORK = 0
       INTEGER :: N, NN, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR
    !  .. LOCAL ARRAYS ..
       REAL(dp), TARGET :: LLVL(1,1), LLVR(1,1)
       REAL(dp), POINTER :: WORK(:)
    !  .. INTRINSIC FUNCTIONS ..
       INTRINSIC MAX, PRESENT, SIZE
    !  .. EXECUTABLE STATEMENTS ..
       LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
       IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
       ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
       IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
       ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
    !  .. TEST THE ARGUMENTS
       IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
       ELSE IF( SIZE( WR ) /= N )THEN; LINFO = -2
       ELSE IF( SIZE( WI ) /= N )THEN; LINFO = -3
       ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -4
       ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -5
       ELSE IF( N > 0 )THEN
          NN = 3; IF( LSAME(LJOBVL,'V').OR.LSAME(LJOBVR,'V') ) NN = NN + 1
          LWORK = MAX( 1, NN*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
             IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
                LWORK = MAX( 1, NN*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
                IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
             END IF
          IF( ISTAT == 0 ) THEN
            IF( PRESENT(VL) )THEN
               IF( PRESENT(VR) )THEN
                   CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            VL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
           ELSE
              CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            VL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
           ENDIF
         ELSE
           IF( PRESENT(VR) )THEN
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            LLVL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
           ELSE
               CALL la_GEEV_m( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                            LLVL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
           ENDIF
         ENDIF  
             IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
          ELSE; LINFO = -100; ENDIF
          DEALLOCATE(WORK, STAT=ISTAT1)
       ENDIF
       CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
    end subroutine std_lapack_geev_real
    
    
    subroutine std_lapack_gesv_real( A, B, IPIV, INFO )
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      REAL(dp), INTENT(INOUT) :: A(:,:), B(:,:)

!     .. "Parameters" ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. "Local Pointers" ..
      INTEGER, POINTER :: LPIV(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; ISTAT = 0
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     .. "Test the arguments" ..
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN
         LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= SIZE(A,1) .OR. SIZE(B,2) < 0 ) THEN
         LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) )THEN
            LINFO = -3
      ELSE IF ( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) )THEN
            LPIV => IPIV
         ELSE
            ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT )
         END IF
         IF ( ISTAT == 0 ) THEN
!        .. "Call LAPACK77 routine" ..
            CALL la_GESV_m( SIZE(A,1), SIZE(B,2), A, MAX(1,SIZE(A,1)), LPIV, B, MAX(1,SIZE(A,1)), &
                           LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(IPIV) )THEN
            DEALLOCATE(LPIV, STAT = ISTAT1 )
         END IF
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )

    end subroutine std_lapack_gesv_real

    subroutine std_lapack_gesv_cmplx( A, B, IPIV, INFO )
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      COMPLEX(dpc), INTENT(INOUT) :: A(:,:), B(:,:)

!     .. "Parameters" ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. "Local Pointers" ..
      INTEGER, POINTER :: LPIV(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; ISTAT = 0
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     .. "Test the arguments" ..
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN
         LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= SIZE(A,1) .OR. SIZE(B,2) < 0 ) THEN
         LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) )THEN
            LINFO = -3
      ELSE IF ( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) )THEN
            LPIV => IPIV
         ELSE
            ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT )
         END IF
         IF ( ISTAT == 0 ) THEN
!        .. "Call LAPACK77 routine" ..
            CALL la_GESV_m( SIZE(A,1), SIZE(B,2), A, MAX(1,SIZE(A,1)), LPIV, B, MAX(1,SIZE(A,1)), &
                           LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(IPIV) )THEN
            DEALLOCATE(LPIV, STAT = ISTAT1 )
         END IF
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )

    end subroutine std_lapack_gesv_cmplx


     LOGICAL FUNCTION LSAME( CA, CB )
!
!  PURPOSE
!  =======
!
!  LSAME  TESTS IF CA IS THE SAME LETTER AS CB REGARDLESS OF CASE.
!
!  PARAMETERS
!  ==========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CHARACTERS TO BE COMPARED.
!
!  .. SCALAR ARGUMENTS ..
      CHARACTER*1, INTENT(IN) :: CA, CB
!  .. PARAMETERS ..
      INTEGER, PARAMETER      :: IOFF=32
!  .. LOCAL SCALARS ..
      INTEGER                 :: INTA, INTB, ZCODE
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC                  ICHAR
!
!  .. EXECUTABLE STATEMENTS ..
!
!  TEST IF THE CHARACTERS ARE EQUAL
!
      LSAME = CA == CB
!
!  NOW TEST FOR EQUIVALENCE
!
      IF( .NOT.LSAME )THEN
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
         ZCODE = ICHAR( 'Z' )
!
         INTA = ICHAR( CA )
         INTB = ICHAR( CB )
!
         IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 )THEN
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
            IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
            IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
         ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 )THEN
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
!    &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
         ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 )THEN
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
            IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
         ENDIF
         LSAME = INTA == INTB
      ENDIF
      END FUNCTION LSAME

      SUBROUTINE ERINFO(LINFO, SRNAME, INFO, ISTAT)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. IMPLICIT STATEMENT ..
         IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
         CHARACTER( LEN = * ), INTENT(IN)              :: SRNAME
         INTEGER             , INTENT(IN)              :: LINFO
         INTEGER             , INTENT(OUT), OPTIONAL   :: INFO
         INTEGER             , INTENT(IN), OPTIONAL    :: ISTAT
!  .. EXECUTABLE STATEMENTS ..
!         IF( ( LINFO < 0 .AND. LINFO > -200 ) .OR.                     &
!    &       ( LINFO > 0 .AND. .NOT.PRESENT(INFO) ) )THEN
      IF( ( ( LINFO < 0 .AND. LINFO > -200 ) .OR. LINFO > 0 )           &
     &           .AND. .NOT.PRESENT(INFO) )THEN
        WRITE (*,*) 'Program terminated in LAPACK95 subroutine ',SRNAME
        WRITE (*,*) 'Error indicator, INFO = ',LINFO
        IF( PRESENT(ISTAT) )THEN
          IF( ISTAT /= 0 ) THEN
            IF( LINFO == -100 )THEN
              WRITE (*,*) 'The statement ALLOCATE causes STATUS = ',    &
     &                    ISTAT
            ELSE
              WRITE (*,*) 'LINFO = ', LINFO, ' not expected'
            END IF
          END IF   
        END IF
        STOP
         ELSE IF( LINFO <= -200 ) THEN
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
           WRITE(*,*) '*** WARNING, INFO = ', LINFO, ' WARNING ***'
           IF( LINFO == -200 )THEN
             WRITE(*,*)                                                 &
     &        'Could not allocate sufficient workspace for the optimum'
             WRITE(*,*)                                                 &
     &        'blocksize, hence the routine may not have performed as'
             WRITE(*,*) 'efficiently as possible'
         ELSE
           WRITE(*,*) 'Unexpected warning'
         END IF
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
        END IF
        IF( PRESENT(INFO) ) THEN
          INFO = LINFO
        END IF
      END SUBROUTINE ERINFO


      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      REAL(DP)   AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPBTRF computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  Contributed by
!  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   WORK( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DPBTF2, DPOTF2, DSYRK, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
         ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
        RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KD ) THEN
!
!        Use unblocked code
!
         CALL DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A12
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose',              &
                                'Non-unit', IB, I2, ONE, AB( KD+1, I ),     &
                                LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Upper', 'Transpose', I2, IB, -ONE,    &
                                 AB( KD+1-IB, I+IB ), LDAB-1, ONE,      &
                                 AB( KD+1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
!
!                    Update A13 (in the work array).
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose',              &
                                 'Non-unit', IB, I3, ONE, AB( KD+1, I ),    &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                     IF( I2.GT.0 )                                          &
                        CALL DGEMM( 'Transpose', 'No Transpose', I2, I3,    &
                                    IB, -ONE, AB( KD+1-IB, I+IB ),          &
                                    LDAB-1, WORK, LDWORK, ONE,              &
                                    AB( 1+IB, I+KD ), LDAB-1 )
!
!                    Update A33
!
                     CALL DSYRK( 'Upper', 'Transpose', I3, IB, -ONE,        &
                                 WORK, LDWORK, ONE, AB( KD+1, I+KD ),       &
                                 LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A21
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose',         &
                                 'Non-unit', I2, IB, ONE, AB( 1, I ),   &
                                 LDAB-1, AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Lower', 'No Transpose', I2, IB, -ONE, &
                                 AB( 1+IB, I ), LDAB-1, ONE,            &
                                 AB( 1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
!
!                    Update A31 (in the work array).
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose',         &
                                 'Non-unit', I3, IB, ONE, AB( 1, I ),   &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A32
!
                     IF( I2.GT.0 )                                      &
                        CALL DGEMM( 'No transpose', 'Transpose', I3, I2,&
                                    IB, -ONE, WORK, LDWORK,             &
                                    AB( 1+IB, I ), LDAB-1, ONE,         &
                                    AB( 1+KD-IB, I+IB ), LDAB-1 )
!
!                    Update A33
!
                     CALL DSYRK( 'Lower', 'No Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( 1, I+KD ),      &
                                 LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      RETURN
!
  150 CONTINUE
      RETURN
!
!     End of DPBTRF
!
      END SUBROUTINE

      SUBROUTINE ZPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
!     ..
!     .. Array Arguments ..
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZPBTRF computes the Cholesky factorization of a complex Hermitian
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U**H * U,  if UPLO = 'U', or
!     A = L  * L**H,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the Hermitian band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**H*U or A = L*L**H of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  Contributed by
!  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, IB, II, J, JJ, NB
!     ..
!     .. Local Arrays ..
      COMPLEX*16         WORK( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZHERK, ZPBTF2, ZPOTF2, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND.  &
          ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
        RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'ZPBTRF', UPLO, N, KD, -1, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KD ) THEN
!
!        Use unblocked code
!
         CALL ZPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
!
!        Use blocked code
!
         IF( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute the Cholesky factorization of a Hermitian band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL ZPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A12
!
                     CALL ZTRSM( 'Left', 'Upper', 'Conjugate transpose',&
                                'Non-unit', IB, I2, CONE,&
                                AB( KD+1, I ), LDAB-1,&
                                AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                     CALL ZHERK( 'Upper', 'Conjugate transpose', I2, IB,&
                                -ONE, AB( KD+1-IB, I+IB ), LDAB-1, ONE,&
                                AB( KD+1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE
!
!                    Update A13 (in the work array).
!
                     CALL ZTRSM( 'Left', 'Upper', 'Conjugate transpose',&
                                'Non-unit', IB, I3, CONE,&
                                AB( KD+1, I ), LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                     IF( I2.GT.0 )&
                       CALL ZGEMM( 'Conjugate transpose',&
                                   'No transpose', I2, I3, IB, -CONE,&
                                   AB( KD+1-IB, I+IB ), LDAB-1, WORK,&
                                   LDWORK, CONE, AB( 1+IB, I+KD ),&
                                   LDAB-1 )
!
!                    Update A33
!
                     CALL ZHERK( 'Upper', 'Conjugate transpose', I3, IB,&
                                -ONE, WORK, LDWORK, ONE,&
                                AB( KD+1, I+KD ), LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
!
!           Compute the Cholesky factorization of a Hermitian band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE
!
!           Process the band matrix one diagonal block at a time.
!
            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL ZPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II.NE.0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB.LE.N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A21
!
                     CALL ZTRSM( 'Right', 'Lower',&
                                'Conjugate transpose', 'Non-unit', I2,&
                                IB, CONE, AB( 1, I ), LDAB-1,&
                                AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                     CALL ZHERK( 'Lower', 'No transpose', I2, IB, -ONE,&
                                AB( 1+IB, I ), LDAB-1, ONE,&
                                AB( 1, I+IB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE
!
!                    Update A31 (in the work array).
!
                     CALL ZTRSM( 'Right', 'Lower',&
                                'Conjugate transpose', 'Non-unit', I3,&
                                IB, CONE, AB( 1, I ), LDAB-1, WORK,&
                                LDWORK )
!
!                    Update A32
!
                     IF( I2.GT.0 )&
                       CALL ZGEMM( 'No transpose',&
                                   'Conjugate transpose', I3, I2, IB,&
                                   -CONE, WORK, LDWORK, AB( 1+IB, I ),&
                                   LDAB-1, CONE, AB( 1+KD-IB, I+IB ),&
                                   LDAB-1 )
!
!                    Update A33
!
                     CALL ZHERK( 'Lower', 'No transpose', I3, IB, -ONE,&
                                WORK, LDWORK, ONE, AB( 1, I+KD ),&
                                LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      RETURN
!
  150 CONTINUE
      RETURN
!
!     End of ZPBTRF
!
      END SUBROUTINE

    
end module std_lapack
