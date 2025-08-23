
#include "panic.F90"

module numa__dlahqr
	implicit none

	! TODO: DRY

	interface outer_product
		!procedure :: outer_product_c64
		procedure :: outer_product_f64
	end interface outer_product

contains

function outer_product_f64(a, b) result(c)
	double precision, intent(in) :: a(:), b(:)
	double precision, allocatable :: c(:,:)

	integer :: i, j, na, nb
	na = size(a)
	nb = size(b)

	allocate(c(na, nb))
	do j = 1, nb
	do i = 1, na
		c(i,j) = a(i) * b(j)
	end do
	end do

end function outer_product_f64

function eye(n)
	! n x n identity matrix

	integer, intent(in) :: n
	double precision, allocatable :: eye(:,:)

	integer :: i, j
	allocate(eye(n, n))
	do i = 1, n
		do j = 1, n
			if (i == j) then
				eye(i,j) = 1.d0
			else
				eye(i,j) = 0.d0
			end if
		end do
	end do

end function eye

!********

double precision function sign_(x)
	double precision, intent(in) :: x
	sign_ = sign(1.d0, x)
end function sign_

!********

function house(x, iostat) result(pp)
	! Return the Householder reflector `pp` such that pp * x == [1, 0, 0, ...]
	!
	! Could just return `v` like house_c64(), but this fn is only used for 3x3
	! matrices, whereas house_c64() runs on arbitrarily large matrices for QR
	! factoring

	use numa__utils
	double precision, intent(in) :: x(:)
	double precision, allocatable :: pp(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision :: alpha, normv
	double precision, allocatable :: v(:)
	integer :: n

	if (present(iostat)) iostat = 0

	n = size(x)

	alpha = -sign_(x(1)) * norm2(x)
	v = x
	v(1) = v(1) - alpha

	normv = norm2(v)
	if (normv <= 0) then
		v = v * 1.d50
		normv = norm2(v)
		print *, "normv = ", normv
		print *, "v/normv = ", v / normv

		pp = eye(n)
		print *, "v = ", v
		msg = "vector is singular in house()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	v = v / normv

	!print *, "v house = ", v

	pp = eye(n) - 2.d0 * outer_product(v, v)

	!print *, "in house():"
	!print *, "x = ", x
	!print *, "pp = "
	!print "(3es15.5)", pp
	!print *, "pp * x = ", matmul(pp, x)

end function house

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
                         SAFMIN, SMLNUM, SN, T1, TR, TST, &
                         ULP
      INTEGER            I, I1, I2, i3, ITS, ITMAX, J, K, L, M, NH, NR, NZ, &
                         KDEFL 
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
	  double precision :: p2(2,2), p3(3,3)
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
      SAFMIN = tiny(safmin)
      SAFMAX = ONE / SAFMIN
      ulp = epsilon(ulp)

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
            !print *, "nr = ", nr
            IF( K.GT.M ) &
               v(1: nr) = h(k: k+nr-1, k-1)

            !print *, "v       = ", v

            !! TODO: delete if unused
            !p3 = house(v)

            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            !print *, "v*tau   = ", v * t1
            !print *, "v       = ", v
            !print *, "t1 = ", t1

            IF( NR.EQ.3 ) THEN
               p3 = eye(nr) - t1 * outer_product([1.d0, v(2:nr)], [1.d0, v(2:nr)])

!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
               h(k:k+2, k:i2) = matmul(p3, h(k:k+2, k:i2))

!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
               i3 = min(k+3, i)
               h(i1:i3, k:k+2) = matmul(h(i1:i3, k:k+2), p3)
!
               IF( WANTZ ) THEN
!                 Accumulate transformations in the matrix Z
                  z(iloz:ihiz, k:k+2) = matmul(z(iloz:ihiz, k:k+2), p3)
               END IF
            ELSE IF( NR.EQ.2 ) THEN
				!print *, "nr == 2"
               p2 = eye(nr) - t1 * outer_product([1.d0, v(2:nr)], [1.d0, v(2:nr)])

!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
               h(k:k+1, k:i2) = matmul(p2, h(k:k+1, k:i2))

!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
               h(i1:i, k:k+1) = matmul(h(i1:i, k:k+1), p2)
               IF( WANTZ ) THEN
!                 Accumulate transformations in the matrix Z
                  z(iloz:ihiz, k:k+1) = matmul(z(iloz:ihiz, k:k+1), p2)
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
         if (wantt .or. wantz) then
            p2(1,:) = [ cs, sn]
            p2(2,:) = [-sn, cs]
         end if

         IF( WANTT ) THEN
!           Apply the transformation to the rest of H.
            IF( I2.GT.I ) then
               h(i-1:i, i+1:i2) = matmul(p2, h(i-1:i, i+1:i2))
            end if
            h(i1: i-2, i-1:i) = matmul(h(i1: i-2, i-1:i), transpose(p2))
         END IF
         IF( WANTZ ) THEN
!           Apply the transformation to Z.
            z(iloz:nz, i-1:i) = matmul(z(iloz:nz, i-1:i), transpose(p2))
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
!
      END SUBROUTINE DLAHQR
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
      SAFMIN = tiny(safmin)
      eps = epsilon(eps)
      SAFMN2 = radix(safmn2)**INT( LOG( SAFMIN / EPS ) / &
                  LOG( dble(radix(safmn2)) ) / TWO )

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
            TAU = norm2([C, Z])
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
            TAU = norm2([SIGMA, TEMP])
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
      END SUBROUTINE DLANV2
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
      XNORM = norm2(x(1: n-1))
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
         BETA = -SIGN(norm2([ALPHA, XNORM]), ALPHA )
         SAFMIN = tiny(safmin) / (epsilon(safmin) / radix(safmin))
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1

            x(1: n-1: incx) = rsafmn * x(1: n-1: incx)
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = norm2(x(1: n-1))
            BETA = -SIGN(norm2([ALPHA, XNORM]), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         x(1: n-1: incx) = x(1: n-1: incx) / (alpha - beta)
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
!
      END SUBROUTINE DLARFG
!***********************************************************************
end module numa__dlahqr
