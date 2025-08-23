
#include "panic.F90"

module numa__dlahqr
	implicit none
contains

!> \brief \b dlahqr computes the eigenvalues and schur factorization of an upper hessenberg matrix, using the double-shift/single-shift qr algorithm.

!  =========== documentation ===========

! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/

!> \htmlonly
!> download dlahqr + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlahqr.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlahqr.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlahqr.f">
!> [txt]</a>
!> \endhtmlonly

!  definition:
!  ===========

!       subroutine dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi,
!                          iloz, ihiz, z, ldz, info )

!       .. scalar arguments ..
!       integer            ihi, ihiz, ilo, iloz, info, ldh, ldz, n
!       logical            wantt, wantz
!       ..
!       .. array arguments ..
!       double precision   h( ldh, * ), wi( * ), wr( * ), z( ldz, * )
!       ..


!> \par purpose:
!  =============
!>
!> \verbatim
!>
!>    dlahqr is an auxiliary routine called by dhseqr to update the
!>    eigenvalues and schur decomposition already computed by dhseqr, by
!>    dealing with the hessenberg submatrix in rows and columns ilo to
!>    ihi.
!> \endverbatim

!  arguments:
!  ==========

!> \param[in] wantt
!> \verbatim
!>          wantt is logical
!>          = .true. : the full schur form t is required;
!>          = .false.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] wantz
!> \verbatim
!>          wantz is logical
!>          = .true. : the matrix of schur vectors z is required;
!>          = .false.: schur vectors are not required.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the matrix h.  n >= 0.
!> \endverbatim
!>
!> \param[in] ilo
!> \verbatim
!>          ilo is integer
!> \endverbatim
!>
!> \param[in] ihi
!> \verbatim
!>          ihi is integer
!>          it is assumed that h is already upper quasi-triangular in
!>          rows and columns ihi+1:n, and that h(ilo,ilo-1) = 0 (unless
!>          ilo = 1). dlahqr works primarily with the hessenberg
!>          submatrix in rows and columns ilo to ihi, but applies
!>          transformations to all of h if wantt is .true..
!>          1 <= ilo <= max(1,ihi); ihi <= n.
!> \endverbatim
!>
!> \param[in,out] h
!> \verbatim
!>          h is double precision array, dimension (ldh,n)
!>          on entry, the upper hessenberg matrix h.
!>          on exit, if info is zero and if wantt is .true., h is upper
!>          quasi-triangular in rows and columns ilo:ihi, with any
!>          2-by-2 diagonal blocks in standard form. if info is zero
!>          and wantt is .false., the contents of h are unspecified on
!>          exit.  the output state of h if info is nonzero is given
!>          below under the description of info.
!> \endverbatim
!>
!> \param[in] ldh
!> \verbatim
!>          ldh is integer
!>          the leading dimension of the array h. ldh >= max(1,n).
!> \endverbatim
!>
!> \param[out] wr
!> \verbatim
!>          wr is double precision array, dimension (n)
!> \endverbatim
!>
!> \param[out] wi
!> \verbatim
!>          wi is double precision array, dimension (n)
!>          the real and imaginary parts, respectively, of the computed
!>          eigenvalues ilo to ihi are stored in the corresponding
!>          elements of wr and wi. if two eigenvalues are computed as a
!>          complex conjugate pair, they are stored in consecutive
!>          elements of wr and wi, say the i-th and (i+1)th, with
!>          wi(i) > 0 and wi(i+1) < 0. if wantt is .true., the
!>          eigenvalues are stored in the same order as on the diagonal
!>          of the schur form returned in h, with wr(i) = h(i,i), and, if
!>          h(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!>          wi(i) = sqrt(h(i+1,i)*h(i,i+1)) and wi(i+1) = -wi(i).
!> \endverbatim
!>
!> \param[in] iloz
!> \verbatim
!>          iloz is integer
!> \endverbatim
!>
!> \param[in] ihiz
!> \verbatim
!>          ihiz is integer
!>          specify the rows of z to which transformations must be
!>          applied if wantz is .true..
!>          1 <= iloz <= ilo; ihi <= ihiz <= n.
!> \endverbatim
!>
!> \param[in,out] z
!> \verbatim
!>          z is double precision array, dimension (ldz,n)
!>          if wantz is .true., on entry z must contain the current
!>          matrix z of transformations accumulated by dhseqr, and on
!>          exit z has been updated; transformations are applied only to
!>          the submatrix z(iloz:ihiz,ilo:ihi).
!>          if wantz is .false., z is not referenced.
!> \endverbatim
!>
!> \param[in] ldz
!> \verbatim
!>          ldz is integer
!>          the leading dimension of the array z. ldz >= max(1,n).
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is integer
!>           = 0:  successful exit
!>           > 0:  if info = i, dlahqr failed to compute all the
!>                  eigenvalues ilo to ihi in a total of 30 iterations
!>                  per eigenvalue; elements i+1:ihi of wr and wi
!>                  contain those eigenvalues which have been
!>                  successfully computed.
!>
!>                  if info > 0 and wantt is .false., then on exit,
!>                  the remaining unconverged eigenvalues are the
!>                  eigenvalues of the upper hessenberg matrix rows
!>                  and columns ilo through info of the final, output
!>                  value of h.
!>
!>                  if info > 0 and wantt is .true., then on exit
!>          (*)       (initial value of h)*u  = u*(final value of h)
!>                  where u is an orthogonal matrix.    the final
!>                  value of h is upper hessenberg and triangular in
!>                  rows and columns info+1 through ihi.
!>
!>                  if info > 0 and wantz is .true., then on exit
!>                      (final value of z)  = (initial value of z)*u
!>                  where u is the orthogonal matrix in (*)
!>                  (regardless of the value of wantt.)
!> \endverbatim

!  authors:
!  ========

!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.

!> \ingroup lahqr

!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>     02-96 based on modifications by
!>     david day, sandia national laboratory, usa
!>
!>     12-04 further modifications by
!>     ralph byers, university of kansas, usa
!>     this is a modified version of dlahqr from lapack version 3.0.
!>     it is (1) more robust against overflow and underflow and
!>     (2) adopts the more conservative ahues & tisseur stopping
!>     criterion (lawn 122, 1997).
!> \endverbatim
!>
!  =====================================================================
      subroutine dlahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, &
                         iloz, ihiz, z, ldz, info )

      use numa__utils
      implicit none

!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--

!     .. scalar arguments ..
      integer            ihi, ihiz, ilo, iloz, info, ldh, ldz, n
      logical            wantt, wantz
!     ..
!     .. array arguments ..
      double precision   h( ldh, * ), wi( * ), wr( * ), z( ldz, * )
!     ..

!  =========================================================

!     .. parameters ..
	  double precision, parameter :: DAT1 = 3.0d0 / 4.0d0, DAT2 = -0.4375d0
      integer, parameter :: KEXSH = 10
!     ..
!     .. local scalars ..
      double precision   aa, ab, ba, bb, cs, h11, h12, h21, h21s, &
                         h22, rt1i, rt1r, rt2i, rt2r, rtdisc, s, safmax, &
                         safmin, smlnum, sn, t1, tr, tst, &
                         ulp, det_
      integer            i, i1, i2, i3, its, itmax, k, l, m, nh, nr, nz, &
                         kdefl
      logical :: has_converged
!     ..
!     .. local arrays ..
      double precision   v( 3 )
	  double precision :: p2(2,2), p3(3,3)
!     ..
!     ..
!     .. intrinsic functions ..
      intrinsic          abs, dble, max, min, sqrt
!     ..
!     .. executable statements ..

      info = 0

!     quick return if possible

      if( n == 0 ) &
         return
      if( ilo == ihi ) then
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = 0
         return
      end if

      if( ilo <= ihi-2 ) &
         h( ihi, ihi-2 ) = 0

      nh = ihi - ilo + 1
      nz = ihiz - iloz + 1

!     set machine-dependent constants for the stopping criterion.
      safmin = tiny(safmin)
      safmax = 1 / safmin
      ulp = epsilon(ulp)

      smlnum = safmin*( dble( nh ) / ulp )

!     i1 and i2 are the indices of the first row and last column of h
!     to which transformations must be applied. if eigenvalues only are
!     being computed, i1 and i2 are set inside the main loop.

      if( wantt ) then
         i1 = 1
         i2 = n
      end if

!     itmax is the total number of qr iterations allowed.

      itmax = 30 * max( 10, nh )

!     kdefl counts the number of iterations since a deflation

      kdefl = 0

!     the main loop begins here. i is the loop index and decreases from
!     ihi to ilo in steps of 1 or 2. each iteration of the loop works
!     with the active submatrix in rows and columns l to i.
!     eigenvalues i+1 to ihi have already converged. either l = ilo or
!     h(l,l-1) is negligible so that the matrix splits.
      i = ihi
      do while (i >= ilo)
      l = ilo

!     perform qr iterations on rows and columns ilo to i until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.

      has_converged = .false.
      do its = 0, itmax

!        look for a single small subdiagonal element.

         do k = i, l + 1, -1
            if( abs( h( k, k-1 ) ) <= smlnum ) &
               exit
            tst = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            if( tst == 0 ) then
               if( k-2 >= ilo ) &
                  tst = tst + abs( h( k-1, k-2 ) )
               if( k+1 <= ihi ) &
                  tst = tst + abs( h( k+1, k ) )
            end if
!           ==== the following is a conservative small subdiagonal
!           .    deflation  criterion due to ahues & tisseur (lawn 122,
!           .    1997). it has better mathematical foundation and
!           .    improves accuracy in some cases.  ====
            if( abs( h( k, k-1 ) ) <= ulp*tst ) then
               ab = max( abs( h( k, k-1 ) ), abs( h( k-1, k ) ) )
               ba = min( abs( h( k, k-1 ) ), abs( h( k-1, k ) ) )
               aa = max( abs( h( k, k ) ), &
                    abs( h( k-1, k-1 )-h( k, k ) ) )
               bb = min( abs( h( k, k ) ), &
                    abs( h( k-1, k-1 )-h( k, k ) ) )
               s = aa + ab
               if( ba*( ab / s ) <= max( smlnum, &
                   ulp*( bb*( aa / s ) ) ) ) exit
            end if
		 end do
         l = k
         if( l > ilo ) then

!           h(l,l-1) is negligible

            h( l, l-1 ) = 0
         end if

!        exit from loop if a submatrix of order 1 or 2 has split off.
         if( l >= i-1 ) then
            has_converged = .true.
            exit
         end if
         kdefl = kdefl + 1

!        now the active submatrix is in rows and columns l to i. if
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
         if( .not.wantt ) then
            i1 = l
            i2 = i
         end if

         if( mod(kdefl,2*KEXSH) == 0 ) then
!           exceptional shift.
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h11 = DAT1*s + h( i, i )
            h12 = DAT2*s
            h21 = s
            h22 = h11
         else if( mod(kdefl,KEXSH) == 0 ) then
!           exceptional shift.
            s = abs( h( l+1, l ) ) + abs( h( l+2, l+1 ) )
            h11 = DAT1*s + h( l, l )
            h12 = DAT2*s
            h21 = s
            h22 = h11
         else
!           prepare to use francis' double shift
!           (i.e. 2nd degree generalized rayleigh quotient)
            h11 = h( i-1, i-1 )
            h21 = h( i, i-1 )
            h12 = h( i-1, i )
            h22 = h( i, i )
         end if
         s = abs( h11 ) + abs( h12 ) + abs( h21 ) + abs( h22 )
         if( s <= 0 ) then
            rt1r = 0
            rt1i = 0
            rt2r = 0
            rt2i = 0
         else
            h11 = h11 / s
            h21 = h21 / s
            h12 = h12 / s
            h22 = h22 / s
            tr = ( h11+h22 ) / 2
            det_ = ( h11-tr )*( h22-tr ) - h12*h21
            rtdisc = sqrt( abs( det_ ) )
            if( det_ >= 0 ) then

!              ==== complex conjugate shifts ====

               rt1r = tr*s
               rt2r = rt1r
               rt1i = rtdisc*s
               rt2i = -rt1i
            else

!              ==== real shifts (use only one of them)  ====

               rt1r = tr + rtdisc
               rt2r = tr - rtdisc
               if( abs( rt1r-h22 ) <= abs( rt2r-h22 ) ) then
                  rt1r = rt1r*s
                  rt2r = rt1r
               else
                  rt2r = rt2r*s
                  rt1r = rt2r
               end if
               rt1i = 0
               rt2i = 0
            end if
         end if

!        look for two consecutive small subdiagonal elements.

         do m = i - 2, l, -1
!           determine the effect of starting the double-shift qr
!           iteration at row m, and see if this would make h(m,m-1)
!           negligible.  (the following uses scaling to avoid
!           overflows and most underflows.)

            h21s = h( m+1, m )
            s = abs( h( m, m )-rt2r ) + abs( rt2i ) + abs( h21s )
            h21s = h( m+1, m ) / s
            v( 1 ) = h21s*h( m, m+1 ) + ( h( m, m )-rt1r )* &
                     ( ( h( m, m )-rt2r ) / s ) - rt1i*( rt2i / s )
            v( 2 ) = h21s*( h( m, m )+h( m+1, m+1 )-rt1r-rt2r )
            v( 3 ) = h21s*h( m+2, m+1 )
            s = abs( v( 1 ) ) + abs( v( 2 ) ) + abs( v( 3 ) )
            v( 1 ) = v( 1 ) / s
            v( 2 ) = v( 2 ) / s
            v( 3 ) = v( 3 ) / s
            if( m == l ) &
               exit
            if( abs( h( m, m-1 ) )*( abs( v( 2 ) )+abs( v( 3 ) ) ) <=  &
                ulp*abs( v( 1 ) )*( abs( h( m-1, m-1 ) )+abs( h( m, &
                m ) )+abs( h( m+1, m+1 ) ) ) ) exit
         end do

!        double-shift qr step
         do k = m, i - 1

!           the first iteration of this loop determines a reflection g
!           from the vector v and applies it from left and right to h,
!           thus creating a nonzero bulge below the subdiagonal.

!           each subsequent iteration determines a reflection g to
!           restore the hessenberg form in the (k-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. nr is the order of g.

            nr = min( 3, i-k+1 )
            !print *, "nr = ", nr
            if( k > m ) &
               v(1: nr) = h(k: k+nr-1, k-1)

            !print *, "v       = ", v

            call dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            !print *, "v*tau   = ", v * t1
            !print *, "v       = ", v
            !print *, "t1 = ", t1

            if( k > m ) then
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = 0
               if( k < i-1 ) &
                  h( k+2, k-1 ) = 0
            else if( m > l ) then
!               ==== use the following instead of
!               .    h( k, k-1 ) = -h( k, k-1 ) to
!               .    avoid a bug when v(2) and v(3)
!               .    underflow. ====
               h( k, k-1 ) = h( k, k-1 )*( 1-t1 )
            end if

            if( nr == 3 ) then
               p3 = eye(nr) - t1 * outer_product([1.d0, v(2:nr)], [1.d0, v(2:nr)])

!              apply g from the left to transform the rows of the matrix
!              in columns k to i2.
               h(k:k+2, k:i2) = matmul(p3, h(k:k+2, k:i2))

!              apply g from the right to transform the columns of the
!              matrix in rows i1 to min(k+3,i).
               i3 = min(k+3, i)
               h(i1:i3, k:k+2) = matmul(h(i1:i3, k:k+2), p3)

               if( wantz ) then
!                 accumulate transformations in the matrix z
                  z(iloz:ihiz, k:k+2) = matmul(z(iloz:ihiz, k:k+2), p3)
               end if
            else if( nr == 2 ) then
				!print *, "nr == 2"
               p2 = eye(nr) - t1 * outer_product([1.d0, v(2:nr)], [1.d0, v(2:nr)])

!              apply g from the left to transform the rows of the matrix
!              in columns k to i2.
               h(k:k+1, k:i2) = matmul(p2, h(k:k+1, k:i2))

!              apply g from the right to transform the columns of the
!              matrix in rows i1 to min(k+3,i).
               h(i1:i, k:k+1) = matmul(h(i1:i, k:k+1), p2)
               if( wantz ) then
!                 accumulate transformations in the matrix z
                  z(iloz:ihiz, k:k+1) = matmul(z(iloz:ihiz, k:k+1), p2)
               end if
            end if
         end do
      end do

!     failure to converge in remaining number of iterations
      if (.not. has_converged) then
         info = i
         return
      end if

      if( l == i ) then

!        h(i,i-1) is negligible: one eigenvalue has converged.

         wr( i ) = h( i, i )
         wi( i ) = 0
      else if( l == i-1 ) then
!        h(i-1,i-2) is negligible: a pair of eigenvalues have converged.

!        transform the 2-by-2 submatrix to standard schur form,
!        and compute and store the eigenvalues.

         call dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ), &
                      h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ), &
                      cs, sn )

         if (wantt .or. wantz) then
            p2(1,:) = [ cs, sn]
            p2(2,:) = [-sn, cs]
         end if

         if( wantt ) then
!           apply the transformation to the rest of h.
            if( i2 > i ) then
               h(i-1:i, i+1:i2) = matmul(p2, h(i-1:i, i+1:i2))
            end if
            h(i1: i-2, i-1:i) = matmul(h(i1: i-2, i-1:i), transpose(p2))
         end if
         if( wantz ) then
!           apply the transformation to z.
            z(iloz:nz, i-1:i) = matmul(z(iloz:nz, i-1:i), transpose(p2))
         end if
      end if
!     reset deflation counter
      kdefl = 0

!     return to start of the main loop with new value of i.
      i = l - 1
      end do

      end subroutine dlahqr
!> \brief \b dlanv2 computes the schur factorization of a real 2-by-2 nonsymmetric matrix in standard form.

!  =========== documentation ===========

! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/

!> \htmlonly
!> download dlanv2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlanv2.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlanv2.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlanv2.f">
!> [txt]</a>
!> \endhtmlonly

!  definition:
!  ===========

!       subroutine dlanv2( a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn )

!       .. scalar arguments ..
!       double precision   a, b, c, cs, d, rt1i, rt1r, rt2i, rt2r, sn
!       ..


!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlanv2 computes the schur factorization of a real 2-by-2 nonsymmetric
!> matrix in standard form:
!>
!>      [ a  b ] = [ cs -sn ] [ aa  bb ] [ cs  sn ]
!>      [ c  d ]   [ sn  cs ] [ cc  dd ] [-sn  cs ]
!>
!> where either
!> 1) cc = 0 so that aa and dd are real eigenvalues of the matrix, or
!> 2) aa = dd and bb*cc < 0, so that aa + or - sqrt(bb*cc) are complex
!> conjugate eigenvalues.
!> \endverbatim

!  arguments:
!  ==========

!> \param[in,out] a
!> \verbatim
!>          a is double precision
!> \endverbatim
!>
!> \param[in,out] b
!> \verbatim
!>          b is double precision
!> \endverbatim
!>
!> \param[in,out] c
!> \verbatim
!>          c is double precision
!> \endverbatim
!>
!> \param[in,out] d
!> \verbatim
!>          d is double precision
!>          on entry, the elements of the input matrix.
!>          on exit, they are overwritten by the elements of the
!>          standardised schur form.
!> \endverbatim
!>
!> \param[out] rt1r
!> \verbatim
!>          rt1r is double precision
!> \endverbatim
!>
!> \param[out] rt1i
!> \verbatim
!>          rt1i is double precision
!> \endverbatim
!>
!> \param[out] rt2r
!> \verbatim
!>          rt2r is double precision
!> \endverbatim
!>
!> \param[out] rt2i
!> \verbatim
!>          rt2i is double precision
!>          the real and imaginary parts of the eigenvalues. if the
!>          eigenvalues are a complex conjugate pair, rt1i > 0.
!> \endverbatim
!>
!> \param[out] cs
!> \verbatim
!>          cs is double precision
!> \endverbatim
!>
!> \param[out] sn
!> \verbatim
!>          sn is double precision
!>          parameters of the rotation matrix.
!> \endverbatim

!  authors:
!  ========

!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.

!> \ingroup lanv2

!> \par further details:
!  =====================
!>
!> \verbatim
!>
!>  modified by v. sima, research institute for informatics, bucharest,
!>  romania, to reduce the risk of cancellation errors,
!>  when computing real eigenvalues, and to ensure, if possible, that
!>  abs(rt1r) >= abs(rt2r).
!> \endverbatim
!>
!  =====================================================================
      subroutine dlanv2( a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn )

      use numa__utils

!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--

!     .. scalar arguments ..
      double precision   a, b, c, cs, d, rt1i, rt1r, rt2i, rt2r, sn
!     ..

!  =====================================================================

!     .. parameters ..
      double precision   multpl
      parameter          ( multpl = 4.0d+0 )
!     ..
!     .. local scalars ..
      double precision   aa, bb, bcmax, bcmis, cc, cs1, dd, eps, p, sab, &
                         sac, scale, sigma, sn1, tau, temp, z, safmin, &
                         safmn2, safmx2
      integer            count
      logical :: keep_scaling
!     ..
!     .. executable statements ..

      safmin = tiny(safmin)
      eps = epsilon(eps)
      safmn2 = radix(safmn2)**int( log( safmin / eps ) / &
                  log( dble(radix(safmn2)) ) / 2 )

      safmx2 = 1 / safmn2
      if( c == 0 ) then
         cs = 1
         sn = 0

      else if( b == 0 ) then

!        swap rows and columns

         cs = 0
         sn = 1
         temp = d
         d = a
         a = temp
         b = -c
         c = 0

      else if( ( a-d ) == 0 .and. sign_( b ) /= sign_( c ) ) &
                then
         cs = 1
         sn = 0

      else

         temp = a - d
         p = 0.5d0*temp
         bcmax = max( abs( b ), abs( c ) )
         bcmis = min( abs( b ), abs( c ) )*sign_( b )*sign_( c )
         scale = max( abs( p ), bcmax )
         z = ( p / scale )*p + ( bcmax / scale )*bcmis

!        if z is of the order of the machine accuracy, postpone the
!        decision on the nature of eigenvalues

         if( z >= multpl*eps ) then

!           real eigenvalues. compute a and d.

            z = p + sign( sqrt( scale )*sqrt( z ), p )
            a = d + z
            d = d - ( bcmax / z )*bcmis

!           compute b and the rotation matrix

            tau = norm2([c, z])
            cs = z / tau
            sn = c / tau
            b = b - c
            c = 0

         else

!           complex eigenvalues, or real (almost) equal eigenvalues.
!           make diagonal elements equal.

            count = 0
            sigma = b + c
            keep_scaling = .true.
            do while (keep_scaling)
               keep_scaling = .false.
               count = count + 1
               scale = max( abs(temp), abs(sigma) )
               if( scale >= safmx2 ) then
                  sigma = sigma * safmn2
                  temp = temp * safmn2
                  if (count  <=  20) &
                     keep_scaling = .true.
               end if
               if( scale <= safmn2 ) then
                  sigma = sigma * safmx2
                  temp = temp * safmx2
                  if (count  <=  20) &
                     keep_scaling = .true.
               end if
            end do
            p = 0.5d0*temp
            tau = norm2([sigma, temp])
            cs = sqrt( 0.5d0*( 1+abs( sigma ) / tau ) )
            sn = -( p / ( tau*cs ) )*sign_( sigma )

!           compute [ aa  bb ] = [ a  b ] [ cs -sn ]
!                   [ cc  dd ]   [ c  d ] [ sn  cs ]

            aa = a*cs + b*sn
            bb = -a*sn + b*cs
            cc = c*cs + d*sn
            dd = -c*sn + d*cs

!           compute [ a  b ] = [ cs  sn ] [ aa  bb ]
!                   [ c  d ]   [-sn  cs ] [ cc  dd ]

!           note: some of the multiplications are wrapped in parentheses to
!                 prevent compilers from using fma instructions. see
!                 https://github.com/reference-lapack/lapack/issues/1031.

            a = aa*cs + cc*sn
            b = ( bb*cs ) + ( dd*sn )
            c = -( aa*sn ) + ( cc*cs )
            d = -bb*sn + dd*cs

            temp = 0.5d0*( a+d )
            a = temp
            d = temp

            if( c /= 0 ) then
               if( b /= 0 ) then
                  if( sign_( b ) == sign_( c ) ) then

!                    real eigenvalues: reduce to upper triangular form

                     sab = sqrt( abs( b ) )
                     sac = sqrt( abs( c ) )
                     p = sign( sab*sac, c )
                     tau = 1 / sqrt( abs( b+c ) )
                     a = temp + p
                     d = temp - p
                     b = b - c
                     c = 0
                     cs1 = sab*tau
                     sn1 = sac*tau
                     temp = cs*cs1 - sn*sn1
                     sn = cs*sn1 + sn*cs1
                     cs = temp
                  end if
               else
                  b = -c
                  c = 0
                  temp = cs
                  cs = -sn
                  sn = temp
               end if
            end if
         end if

      end if

!     store eigenvalues in (rt1r,rt1i) and (rt2r,rt2i).

      rt1r = a
      rt2r = d
      if( c == 0 ) then
         rt1i = 0
         rt2i = 0
      else
         rt1i = sqrt( abs( b ) )*sqrt( abs( c ) )
         rt2i = -rt1i
      end if
      end subroutine dlanv2
!> \brief \b dlarfg generates an elementary reflector (householder matrix).

!  =========== documentation ===========

! online html documentation available at
!            http://www.netlib.org/lapack/explore-html/

!> \htmlonly
!> download dlarfg + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f">
!> [tgz]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f">
!> [zip]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f">
!> [txt]</a>
!> \endhtmlonly

!  definition:
!  ===========

!       subroutine dlarfg( n, alpha, x, incx, tau )

!       .. scalar arguments ..
!       integer            incx, n
!       double precision   alpha, tau
!       ..
!       .. array arguments ..
!       double precision   x( * )
!       ..


!> \par purpose:
!  =============
!>
!> \verbatim
!>
!> dlarfg generates a real elementary reflector h of order n, such
!> that
!>
!>       h * ( alpha ) = ( beta ),   h**t * h = i.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. h is represented in the form
!>
!>       h = i - tau * ( 1 ) * ( 1 v**t ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> if the elements of x are all zero, then tau = 0 and h is taken to be
!> the unit matrix.
!>
!> otherwise  1 <= tau <= 2.
!> \endverbatim

!  arguments:
!  ==========

!> \param[in] n
!> \verbatim
!>          n is integer
!>          the order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] alpha
!> \verbatim
!>          alpha is double precision
!>          on entry, the value alpha.
!>          on exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] x
!> \verbatim
!>          x is double precision array, dimension
!>                         (1+(n-2)*abs(incx))
!>          on entry, the vector x.
!>          on exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is integer
!>          the increment between elements of x. incx > 0.
!> \endverbatim
!>
!> \param[out] tau
!> \verbatim
!>          tau is double precision
!>          the value tau.
!> \endverbatim

!  authors:
!  ========

!> \author univ. of tennessee
!> \author univ. of california berkeley
!> \author univ. of colorado denver
!> \author nag ltd.

!> \ingroup larfg

!  =====================================================================
      subroutine dlarfg( n, alpha, x, incx, tau )

!  -- lapack auxiliary routine --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--

!     .. scalar arguments ..
      integer            incx, n
      double precision   alpha, tau
!     ..
!     .. array arguments ..
      double precision   x( * )
!     ..

!  =====================================================================

!     .. local scalars ..
      integer            j, knt
      double precision   beta, rsafmn, safmin, xnorm
!     ..
!     .. executable statements ..

      if( n <= 1 ) then
         tau = 0
         return
      end if

      xnorm = norm2(x(1: n-1))

      if( xnorm == 0 ) then

!        h  =  i

         tau = 0
      else

!        general case

         beta = -sign(norm2([alpha, xnorm]), alpha )
         safmin = tiny(safmin) / (epsilon(safmin) / radix(safmin))
         knt = 0
         if( abs( beta ) < safmin ) then

!           xnorm, beta may be inaccurate; scale x and recompute them
            rsafmn = 1 / safmin
            do
               knt = knt + 1

               x(1: n-1: incx) = rsafmn * x(1: n-1: incx)
               beta = beta*rsafmn
               alpha = alpha*rsafmn
               if( .not. ((abs( beta ) < safmin) .and. (knt  <  20) )) exit
            end do

!           new beta is at most 1, at least safmin
            xnorm = norm2(x(1: n-1))
            beta = -sign(norm2([alpha, xnorm]), alpha )
         end if
         tau = ( beta-alpha ) / beta
         x(1: n-1: incx) = x(1: n-1: incx) / (alpha - beta)

!        if alpha is subnormal, it may lose relative accuracy
         do j = 1, knt
            beta = beta*safmin
         end do
         alpha = beta
      end if

      end subroutine dlarfg
!***********************************************************************
end module numa__dlahqr
