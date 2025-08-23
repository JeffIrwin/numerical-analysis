module numa__dlahqr
contains
!> \brief \b DCOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCOPY copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup copy
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
!
!     End of DCOPY
!
      END
!> \brief \b DISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is DOUBLE PRECISION
!>          Input to test for NaN.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup isnan
!
!  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN
!     ..
!
!  =====================================================================
!
!  ..
!  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
!> \brief \b DLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAHQR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAHQR is an auxiliary routine called by DHSEQR to update the
!>    eigenvalues and Schur decomposition already computed by DHSEQR, by
!>    dealing with the Hessenberg submatrix in rows and columns ILO to
!>    IHI.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          It is assumed that H is already upper quasi-triangular in
!>          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!>          ILO = 1). DLAHQR works primarily with the Hessenberg
!>          submatrix in rows and columns ILO to IHI, but applies
!>          transformations to all of H if WANTT is .TRUE..
!>          1 <= ILO <= max(1,IHI); IHI <= N.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On entry, the upper Hessenberg matrix H.
!>          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
!>          quasi-triangular in rows and columns ILO:IHI, with any
!>          2-by-2 diagonal blocks in standard form. If INFO is zero
!>          and WANTT is .FALSE., the contents of H are unspecified on
!>          exit.  The output state of H if INFO is nonzero is given
!>          below under the description of INFO.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>          The real and imaginary parts, respectively, of the computed
!>          eigenvalues ILO to IHI are stored in the corresponding
!>          elements of WR and WI. If two eigenvalues are computed as a
!>          complex conjugate pair, they are stored in consecutive
!>          elements of WR and WI, say the i-th and (i+1)th, with
!>          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!>          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!>          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE..
!>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          If WANTZ is .TRUE., on entry Z must contain the current
!>          matrix Z of transformations accumulated by DHSEQR, and on
!>          exit Z has been updated; transformations are applied only to
!>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!>          If WANTZ is .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>           = 0:  successful exit
!>           > 0:  If INFO = i, DLAHQR failed to compute all the
!>                  eigenvalues ILO to IHI in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of WR and WI
!>                  contain those eigenvalues which have been
!>                  successfully computed.
!>
!>                  If INFO > 0 and WANTT is .FALSE., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper Hessenberg matrix rows
!>                  and columns ILO through INFO of the final, output
!>                  value of H.
!>
!>                  If INFO > 0 and WANTT is .TRUE., then on exit
!>          (*)       (initial value of H)*U  = U*(final value of H)
!>                  where U is an orthogonal matrix.    The final
!>                  value of H is upper Hessenberg and triangular in
!>                  rows and columns INFO+1 through IHI.
!>
!>                  If INFO > 0 and WANTZ is .TRUE., then on exit
!>                      (final value of Z)  = (initial value of Z)*U
!>                  where U is the orthogonal matrix in (*)
!>                  (regardless of the value of WANTT.)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lahqr
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     02-96 Based on modifications by
!>     David Day, Sandia National Laboratory, USA
!>
!>     12-04 Further modifications by
!>     Ralph Byers, University of Kansas, USA
!>     This is a modified version of DLAHQR from LAPACK version 3.0.
!>     It is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative Ahues & Tisseur stopping
!>     criterion (LAWN 122, 1997).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, INFO )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!     ..
!
!  =========================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0 )
      DOUBLE PRECISION   DAT1, DAT2
      PARAMETER          ( DAT1 = 3.0d0 / 4.0d0, DAT2 = -0.4375d0 )
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 10 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S, &
                         H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX, &
                         SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST, &
                         ULP, V2, V3
      INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ, &
                         KDEFL 
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
!     ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
!
!     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 ) &
         H( IHI, IHI-2 ) = ZERO
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     ITMAX is the total number of QR iterations allowed.
!
      ITMAX = 30 * MAX( 10, NH )
!
!     KDEFL counts the number of iterations since a deflation
!
      KDEFL = 0
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   20 CONTINUE
      L = ILO
      IF( I.LT.ILO ) &
         GO TO 160
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      DO 140 ITS = 0, ITMAX
!
!        Look for a single small subdiagonal element.
!
         DO 30 K = I, L + 1, -1
            IF( ABS( H( K, K-1 ) ).LE.SMLNUM ) &
               GO TO 40
            TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO ) &
                  TST = TST + ABS( H( K-1, K-2 ) )
               IF( K+1.LE.IHI ) &
                  TST = TST + ABS( H( K+1, K ) )
            END IF
!           ==== The following is a conservative small subdiagonal
!           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
!           .    1997). It has better mathematical foundation and
!           .    improves accuracy in some cases.  ====
            IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
               AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               AA = MAX( ABS( H( K, K ) ), &
                    ABS( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( ABS( H( K, K ) ), &
                    ABS( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM, &
                   ULP*( BB*( AA / S ) ) ) )GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
!
!           H(L,L-1) is negligible
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order 1 or 2 has split off.
!
         IF( L.GE.I-1 ) &
            GO TO 150
         KDEFL = KDEFL + 1
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
!
!           Exceptional shift.
!
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H11 = DAT1*S + H( I, I )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
!
!           Exceptional shift.
!
            S = ABS( H( L+1, L ) ) + ABS( H( L+2, L+1 ) )
            H11 = DAT1*S + H( L, L )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE
!
!           Prepare to use Francis' double shift
!           (i.e. 2nd degree generalized Rayleigh quotient)
!
            H11 = H( I-1, I-1 )
            H21 = H( I, I-1 )
            H12 = H( I-1, I )
            H22 = H( I, I )
         END IF
         S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
         IF( S.EQ.ZERO ) THEN
            RT1R = ZERO
            RT1I = ZERO
            RT2R = ZERO
            RT2I = ZERO
         ELSE
            H11 = H11 / S
            H21 = H21 / S
            H12 = H12 / S
            H22 = H22 / S
            TR = ( H11+H22 ) / TWO
            DET = ( H11-TR )*( H22-TR ) - H12*H21
            RTDISC = SQRT( ABS( DET ) )
            IF( DET.GE.ZERO ) THEN
!
!              ==== complex conjugate shifts ====
!
               RT1R = TR*S
               RT2R = RT1R
               RT1I = RTDISC*S
               RT2I = -RT1I
            ELSE
!
!              ==== real shifts (use only one of them)  ====
!
               RT1R = TR + RTDISC
               RT2R = TR - RTDISC
               IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
                  RT1R = RT1R*S
                  RT2R = RT1R
               ELSE
                  RT2R = RT2R*S
                  RT1R = RT2R
               END IF
               RT1I = ZERO
               RT2I = ZERO
            END IF
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 50 M = I - 2, L, -1
!           Determine the effect of starting the double-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.  (The following uses scaling to avoid
!           overflows and most underflows.)
!
            H21S = H( M+1, M )
            S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
            H21S = H( M+1, M ) / S
            V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )* &
                     ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
            V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
            V( 3 ) = H21S*H( M+2, M+1 )
            S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
            V( 1 ) = V( 1 ) / S
            V( 2 ) = V( 2 ) / S
            V( 3 ) = V( 3 ) / S
            IF( M.EQ.L ) &
               GO TO 60
            IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE. &
                ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M, &
                M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
   50    CONTINUE
   60    CONTINUE
!
!        Double-shift QR step
!
         DO 130 K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M ) &
               CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 ) &
                  H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
!               ==== Use the following instead of
!               .    H( K, K-1 ) = -H( K, K-1 ) to
!               .    avoid a bug when v(2) and v(3)
!               .    underflow. ====
               H( K, K-1 ) = H( K, K-1 )*( ONE-T1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 70 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   70          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 80 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   80          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 90 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   90             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 100 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
  100          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 110 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  110          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 120 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
!
  140 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  150 CONTINUE
!
      IF( L.EQ.I ) THEN
!
!        H(I,I-1) is negligible: one eigenvalue has converged.
!
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
!
!        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
!
!        Transform the 2-by-2 submatrix to standard Schur form,
!        and compute and store the eigenvalues.
!
         CALL DLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ), &
                      H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), &
                      CS, SN )
!
         IF( WANTT ) THEN
!
!           Apply the transformation to the rest of H.
!
            IF( I2.GT.I ) &
               CALL DROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH, &
                          CS, SN )
            CALL DROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, &
                       SN )
         END IF
         IF( WANTZ ) THEN
!
!           Apply the transformation to Z.
!
            CALL DROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, &
                       SN )
         END IF
      END IF
!     reset deflation counter
      KDEFL = 0
!
!     return to start of the main loop with new value of I.
!
      I = L - 1
      GO TO 20
!
  160 CONTINUE
      RETURN
!
!     End of DLAHQR
!
      END
!> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is not for general use.  It exists solely to avoid
!> over-optimization in DISNAN.
!>
!> DLAISNAN checks for NaNs by comparing its two arguments for
!> inequality.  NaN is the only floating-point value where NaN != NaN
!> returns .TRUE.  To check for NaNs, pass the same variable as both
!> arguments.
!>
!> A compiler must assume that the two arguments are
!> not the same variable, and the test will not be optimized away.
!> Interprocedural or whole-program optimization may delete this
!> test.  The ISNAN functions will be replaced by the correct
!> Fortran 03 intrinsic once the intrinsic is widely available.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN1
!> \verbatim
!>          DIN1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DIN2
!> \verbatim
!>          DIN2 is DOUBLE PRECISION
!>          Two numbers to compare for inequality.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup laisnan
!
!  =====================================================================
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
!     ..
!
!  =====================================================================
!
!  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
!> \brief \b DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric matrix in standard form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLANV2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanv2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanv2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanv2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
!> matrix in standard form:
!>
!>      [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
!>      [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
!>
!> where either
!> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
!> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
!> conjugate eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION
!>          On entry, the elements of the input matrix.
!>          On exit, they are overwritten by the elements of the
!>          standardised Schur form.
!> \endverbatim
!>
!> \param[out] RT1R
!> \verbatim
!>          RT1R is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT1I
!> \verbatim
!>          RT1I is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT2R
!> \verbatim
!>          RT2R is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] RT2I
!> \verbatim
!>          RT2I is DOUBLE PRECISION
!>          The real and imaginary parts of the eigenvalues. If the
!>          eigenvalues are a complex conjugate pair, RT1I > 0.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is DOUBLE PRECISION
!>          Parameters of the rotation matrix.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lanv2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified by V. Sima, Research Institute for Informatics, Bucharest,
!>  Romania, to reduce the risk of cancellation errors,
!>  when computing real eigenvalues, and to ensure, if possible, that
!>  abs(RT1R) >= abs(RT2R).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0, &
                           TWO = 2.0D0 )
      DOUBLE PRECISION   MULTPL
      PARAMETER          ( MULTPL = 4.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB, &
                         SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z, SAFMIN, &
                         SAFMN2, SAFMX2
      INTEGER            COUNT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'P' )
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                  LOG( DLAMCH( 'B' ) ) / TWO )
      SAFMX2 = ONE / SAFMN2
      IF( C.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
!
      ELSE IF( B.EQ.ZERO ) THEN
!
!        Swap rows and columns
!
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
!
      ELSE IF( ( A-D ).EQ.ZERO .AND. SIGN( ONE, B ).NE.SIGN( ONE, C ) ) &
                THEN
         CS = ONE
         SN = ZERO
!
      ELSE
!
         TEMP = A - D
         P = HALF*TEMP
         BCMAX = MAX( ABS( B ), ABS( C ) )
         BCMIS = MIN( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C )
         SCALE = MAX( ABS( P ), BCMAX )
         Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS
!
!        If Z is of the order of the machine accuracy, postpone the
!        decision on the nature of eigenvalues
!
         IF( Z.GE.MULTPL*EPS ) THEN
!
!           Real eigenvalues. Compute A and D.
!
            Z = P + SIGN( SQRT( SCALE )*SQRT( Z ), P )
            A = D + Z
            D = D - ( BCMAX / Z )*BCMIS
!
!           Compute B and the rotation matrix
!
            TAU = DLAPY2( C, Z )
            CS = Z / TAU
            SN = C / TAU
            B = B - C
            C = ZERO
!
         ELSE
!
!           Complex eigenvalues, or real (almost) equal eigenvalues.
!           Make diagonal elements equal.
!
            COUNT = 0
            SIGMA = B + C
   10       CONTINUE
            COUNT = COUNT + 1
            SCALE = MAX( ABS(TEMP), ABS(SIGMA) )
            IF( SCALE.GE.SAFMX2 ) THEN
               SIGMA = SIGMA * SAFMN2
               TEMP = TEMP * SAFMN2
               IF (COUNT .LE. 20) &
                  GOTO 10
            END IF
            IF( SCALE.LE.SAFMN2 ) THEN
               SIGMA = SIGMA * SAFMX2
               TEMP = TEMP * SAFMX2
               IF (COUNT .LE. 20) &
                  GOTO 10
            END IF
            P = HALF*TEMP
            TAU = DLAPY2( SIGMA, TEMP )
            CS = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
            SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA )
!
!           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
!                   [ CC  DD ]   [ C  D ] [ SN  CS ]
!
            AA = A*CS + B*SN
            BB = -A*SN + B*CS
            CC = C*CS + D*SN
            DD = -C*SN + D*CS
!
!           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
!                   [ C  D ]   [-SN  CS ] [ CC  DD ]
!
!           Note: Some of the multiplications are wrapped in parentheses to
!                 prevent compilers from using FMA instructions. See
!                 https://github.com/Reference-LAPACK/lapack/issues/1031.
!
            A = AA*CS + CC*SN
            B = ( BB*CS ) + ( DD*SN )
            C = -( AA*SN ) + ( CC*CS )
            D = -BB*SN + DD*CS
!
            TEMP = HALF*( A+D )
            A = TEMP
            D = TEMP
!
            IF( C.NE.ZERO ) THEN
               IF( B.NE.ZERO ) THEN
                  IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
!
!                    Real eigenvalues: reduce to upper triangular form
!
                     SAB = SQRT( ABS( B ) )
                     SAC = SQRT( ABS( C ) )
                     P = SIGN( SAB*SAC, C )
                     TAU = ONE / SQRT( ABS( B+C ) )
                     A = TEMP + P
                     D = TEMP - P
                     B = B - C
                     C = ZERO
                     CS1 = SAB*TAU
                     SN1 = SAC*TAU
                     TEMP = CS*CS1 - SN*SN1
                     SN = CS*SN1 + SN*CS1
                     CS = TEMP
                  END IF
               ELSE
                  B = -C
                  C = ZERO
                  TEMP = CS
                  CS = -SN
                  SN = TEMP
               END IF
            END IF
         END IF
!
      END IF
!
!     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
!
      RT1R = A
      RT2R = D
      IF( C.EQ.ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN
!
!     End of DLANV2
!
      END
!> \brief \b DLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAPY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow and unnecessary underflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION
!>          X and Y specify the values x and y.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lapy2
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z, HUGEVAL
      LOGICAL            X_IS_NAN, Y_IS_NAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      X_IS_NAN = DISNAN( X )
      Y_IS_NAN = DISNAN( Y )
      IF ( X_IS_NAN ) DLAPY2 = X
      IF ( Y_IS_NAN ) DLAPY2 = Y
      HUGEVAL = DLAMCH( 'Overflow' )
!
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
            DLAPY2 = W
         ELSE
            DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
!
!     End of DLAPY2
!
      END
!> \brief \b DLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARFG generates a real elementary reflector H of order n, such
!> that
!>
!>       H * ( alpha ) = ( beta ),   H**T * H = I.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> If the elements of x are all zero, then tau = 0 and H is taken to be
!> the unit matrix.
!>
!> Otherwise  1 <= tau <= 2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup larfg
!
!  =====================================================================
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
!
      RETURN
!
!     End of DLARFG
!
      END
!> \brief \b DNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    DNRM2 := sqrt( x'*x )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup nrm2
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
function DNRM2( n, x, incx ) 
   integer, parameter :: wp = kind(1.d0)
   real(wp) :: DNRM2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: maxN = huge(0.0_wp)
!  ..
!  .. Blue's scaling constants ..
   real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
       (minexponent(0._wp) - 1) * 0.5_wp)
   real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
       (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
   real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
       (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
   real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
       (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
   DNRM2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   DNRM2 = scl*sqrt( sumsq )
   return
end function
!> \brief \b DROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION C,S
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DROT applies a plane rotation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup rot
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION C,S
      INTEGER INCX,INCY,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY
!     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!       code for both increments equal to 1
!
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
!
!     End of DROT
!
      END
!> \brief \b DSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DSCAL scales a vector by a constant.
!>    uses unrolled loops for increment equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup scal
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N.LE.0 .OR. INCX.LE.0 .OR. DA.EQ.ONE) RETURN
      IF (INCX.EQ.1) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
!
!        code for increment not equal to 1
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
!
!     End of DSCAL
!
      END
!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION LSAME(CA,CB)
!
!       .. Scalar Arguments ..
!       CHARACTER CA,CB
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!>          CA is CHARACTER*1
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is CHARACTER*1
!>          CA and CB specify the single characters to be compared.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \ingroup lsame
!
!  =====================================================================
      LOGICAL FUNCTION LSAME(CA,CB)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!     .. Scalar Arguments ..
      CHARACTER CA,CB
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
!     ..
!
!     Test if the characters are equal
!
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
!
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
!
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR. &
              INTA.GE.145 .AND. INTA.LE.153 .OR. &
              INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR. &
              INTB.GE.145 .AND. INTB.LE.153 .OR. &
              INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
!
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
!
!     RETURN
!
!     End of LSAME
!
      END
!> \brief \b DLAMCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!     .. Scalar Arguments ..
!     CHARACTER          CMACH
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMCH determines double precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          CMACH is CHARACTER*1
!>          Specifies the value to be returned by DLAMCH:
!>          = 'E' or 'e',   DLAMCH := eps
!>          = 'S' or 's ,   DLAMCH := sfmin
!>          = 'B' or 'b',   DLAMCH := base
!>          = 'P' or 'p',   DLAMCH := eps*base
!>          = 'N' or 'n',   DLAMCH := t
!>          = 'R' or 'r',   DLAMCH := rnd
!>          = 'M' or 'm',   DLAMCH := emin
!>          = 'U' or 'u',   DLAMCH := rmin
!>          = 'L' or 'l',   DLAMCH := emax
!>          = 'O' or 'o',   DLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
!     ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
                         MINEXPONENT, RADIX, TINY
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
      RND = ONE
!
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
!
      DLAMCH = RMACH
      RETURN
!
!     End of DLAMCH
!
      END
!***********************************************************************
!> \brief \b DLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date December 2016
!> \ingroup auxOTHERauxiliary
!>
!> \param[in] A
!> \verbatim
!>          A is a DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is a DOUBLE PRECISION
!>          The values A and B.
!> \endverbatim
!>
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END
!
!***********************************************************************
end module numa__dlahqr
