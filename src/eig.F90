

#include "panic.F90"

!> Module for eigenvalue and eigenvector calculation
module numa__eig

	use numa__blarg
	use numa__linalg

	implicit none

contains

!===============================================================================

subroutine hess(a, u)
	! Apply the Householder Hessenberg reduction to `a`
	!
	! Source:  https://dspace.mit.edu/bitstream/handle/1721.1/75282/18-335j-fall-2006/contents/lecture-notes/lec14.pdf
	use numa__blarg
	use numa__utils
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable, optional, intent(out) :: u(:,:)
	!********

	integer :: i, j, k, n
	double precision :: va
	double precision, allocatable :: v(:), vv(:,:)

	n = size(a,1)
	if (present(u)) allocate(vv(n, n-2))

	do k = 1, n-2

		v = a(k+1:, k)
		v(1) = v(1) + sign_(v(1)) * norm2(v)
		v = v / norm2(v)

		! Avoiding `v` temp array might be a little more work than worthwhile,
		! impossible if `u` is present
		do j = k, n
			va = 2 * dot_product(v, a(k+1:, j))
			a(k+1:, j) = a(k+1:, j) - va * v
		end do
		do i = 1, n
			va = 2 * dot_product(v, a(i, k+1:))
			a(i, k+1:) = a(i, k+1:) - va * v
		end do

		if (present(u)) vv(1: size(v), k) = v

	end do
	if (.not. present(u)) return

	u = eye(n)
	do k = n-2, 1, -1

		v = vv(1: n-k, k)
		do j = k+1, n
			va = 2 * dot_product(v, u(k+1:, j))
			u(k+1:, j) = u(k+1:, j) - va * v
		end do

	end do
	!print *, "u = "
	!print "(4es15.5)", u

end subroutine hess

!===============================================================================

function eig_basic_qr(a, iters, iostat) result(eigvals)
	! Get the real eigenvalues of `a` using `iters` iterations of the basic QR
	! algorithm.  Many iterations may be required for good convergence with
	! large sized `a`
	!
	! The matrix `a` is modified in the process by QR decomposition
	!
	! This only supports real eigenvalues.  The same algorithm could be adapted
	! for complex eigenvalues:
	!
	!     https://math.stackexchange.com/a/3072781/232771
	!
	! The basic idea is that `a` converges to a "quasi-triangular" real matrix
	! in this algorithm.  Where the off-triangle approaches 0, we have real
	! eigenvalues.  Where there is a non-zero in the off-diagonal, we have a
	! complex conjugate pair of eigenvalues, which could be computed as the same
	! eigenvalue of its 2x2 sub-matrix block
	!
	! A better algorithm would also take a tolerance rather than a fixed number
	! of iterations
	!
	! But it's not a good idea to waste time on such improvements for as slow an
	! algorithm as basic QR.  Rather, wait and do the good work on something
	! better like Hessenberg QR, or better yet, a shifting algorithm

	use numa__utils
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable :: eigvals(:)
	integer, intent(in) :: iters
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: diag_(:), q(:,:)
	integer :: i, io

	if (present(iostat)) iostat = 0

	do i = 1, iters
		call qr_factor(a, diag_, iostat = io)
		if (io /= 0) then
			msg = "qr_factor() failed in eig_basic_qr()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		! Could also do a = transpose(qr_mul_transpose(), transpose(r)),
		! but two transposes seems expensive
		q = qr_get_q_expl(a, diag_)
		a = matmul_triu_ge(a, q)
	end do
	!print *, "a = "
	!print "(4es15.5)", a

	eigvals = sorted(diag(a))
	!print *, "eigvals = ", eigvals

end function eig_basic_qr

!===============================================================================

function eig_hess_qr(a, iters, eigvecs, iostat) result(eigvals)
	! Get the real eigenvalues of `a` using `iters` iterations of the Hessenberg
	! QR algorithm

	use numa__blarg
	use numa__utils
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable :: eigvals(:)
	integer, intent(in) :: iters
	double precision, optional, allocatable, intent(out) :: eigvecs(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision :: h1, h2, rad, givens(2,2)
	double precision, allocatable :: c(:), s(:), &
		pq(:,:), r(:,:)
	integer :: i, k, n, io

	if (present(iostat)) iostat = 0
	n = size(a, 1)
	allocate(c(n-1), s(n-1))

	! The matrix `pq` is the product of Q from each step.  Unlike basic QR, we
	! can't initialize pq to eye here.  Instead, we need an initial
	! transformation from the Hessenberg reduction
	call hess(a, pq)

	do i = 1, iters

		do k = 1, n-1
			h1 = a(k, k)
			h2 = a(k+1, k)

			rad = norm2([h1, h2])
			!print *, "rad = ", rad
			c(k) =  h1 / rad
			s(k) = -h2 / rad
			givens(1,:) = [c(k), -s(k)]
			givens(2,:) = [s(k),  c(k)]

			! Ref:  https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
			a(k:k+1, k:) = matmul(givens, a(k:k+1, k:))

		end do
		!print *, "r = "
		!print "(4es15.5)", a

		!! If you stop here, `a` has been overwritten with `r` from its QR
		!! factorization.  You could also compute Q by applying the matmuls in the
		!! loop above to an initial identity matrix, but we don't need Q explicitly
		!! for a Hessenberg QR step
		!stop

		! The rest of the Hessenberg QR step is not part of the QR factorization,
		! rather, it overwrites `a` with R * Q

		! Apply the Givens rotations from the right
		do k = 1, n-1
			givens(1,:) = [ c(k), s(k)]  ! note this is transposed compared to above
			givens(2,:) = [-s(k), c(k)]
			a(1: k+1, k: k+1) = matmul(a(1: k+1, k: k+1), givens)

			if (.not. present(eigvecs)) cycle

			! Update Q product
			pq(:, k:k+1) = matmul(pq(:, k:k+1), givens)

		end do

	end do
	!print *, "a = "
	!print "(4es15.5)", a

	! Don't sort the eigvals, because that will break the ordering of eigvecs
	!eigvals = sorted(diag(a))

	eigvals = diag(a)

	if (.not. present(eigvecs)) return
	!********

	r = triu(a)
	!print *, "r = "
	!print "(4es15.5)", r

	! In order to get the eigenvectors of `a`, we first have to get the
	! eigenvectors of `r` and then transform them by `pq`
	!
	! Ref:  https://math.stackexchange.com/a/3947396/232771
	!
	! Especially the linked code:  https://gist.github.com/uranix/2b4bb821a0e3ffc4531bec547ea67727

	! Find the eigenvectors of `r`
	eigvecs = eye(n)
	do i = 2, n
		! scalar * eye could be optimized here to skip zeros
		eigvecs(1: i-1, i) = invmul(r(i,i) * eye(i-1) - r(:i-1, :i-1), r(:i-1, i), io)

		if (io /= 0) then
			msg = "invmul() failed in eig_hess_qr()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if
	end do
	!print *, "R eigvecs = "
	!print "(4es15.5)", eigvecs

	eigvecs = matmul(pq, eigvecs)
	!print *, "eigvecs hess qr = "
	!print "(4es15.5)", eigvecs

	!print *, "a * w / w ="
	!print "(4es15.5)", matmul(a0, eigvecs) / eigvecs

end function eig_hess_qr

!===============================================================================

function house(x, iostat) result(pp)
	! Return the Householder reflector `pp` such that pp * x == [1, 0, 0, ...]
	!
	! Could just return `v` like house_c64(), but this fn is only used for 3x3
	! matrices, whereas house_c64() runs on arbitrarily large matrices for QR
	! factoring

	use numa__blarg
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
		!print *, "normv = ", normv
		!print *, "v/normv = ", v / normv

		pp = eye(n)
		!print *, "v = ", v
		msg = "vector is singular in house()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	v = v / normv

	pp = eye(n) - 2.d0 * outer_product(v, v)

	!print *, "in house():"
	!print *, "x = ", x
	!print *, "pp = "
	!print "(3es15.5)", pp
	!print *, "pp * x = ", matmul(pp, x)

end function house

!===============================================================================

function eig_lapack(aa, eigvecs, iostat) result(eigvals)
	! Get the eigenvalues of `aa` using the Francis double step QR algorithm,
	! using LAPACK's dlahqr() routine (but still my home-made hess() routine)

	use numa__blarg
	use numa__dlahqr
	use numa__utils
	double precision, intent(inout) :: aa(:,:)
	double complex, allocatable :: eigvals(:)
	double complex, optional, allocatable, intent(out) :: eigvecs(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double complex, allocatable :: ca(:,:), cq(:,:)
	double precision, parameter :: eps = 1.d-10
	double precision, allocatable :: pq(:,:)
	integer :: i, n, io

	if (present(iostat)) iostat = 0

	n = size(aa, 1)

	if (present(eigvecs)) then
		! The matrix `pq` is the product of Q from each step.  Unlike basic QR, we
		! can't initialize pq to eye here.  Instead, we need an initial
		! transformation from the Hessenberg reduction
		call hess(aa, pq)
	else
		call hess(aa)
	end if
	!print *, "aa = "
	!print "("//to_str(n)//"es19.9)", aa

	block
		double precision, allocatable :: wr(:), wi(:)
		integer :: ilo, ihi, ldh, iloz, ihiz, ldz, info
		logical :: wantt, wantz

		! Set indices to factorize all of `aa`, not just a slice
		ilo = 1
		ihi = n
		ldh = n
		iloz = 1
		ihiz = n
		ldz = n

		! Eigval real/imag components
		allocate(wr(n), wi(n))

		if (present(eigvecs)) then
			wantt = .true.  ! want full schur form, not just eigvals
			wantz = .true.
		else
			! Intel requires allocation of `pq` for dlahqr().  For GNU it
			! doesn't matter
			allocate(pq(0,0))
			wantt = .false.
			wantz = .false.
		end if

		call dlahqr(wantt, wantz, n, ilo, ihi, aa, ldh, wr, wi, &
			iloz, ihiz, pq, ldz, info)
		if (info /= 0) then
			msg = "dlahqr() failed in eig_lapack()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		! Copy eigvals out
		eigvals = dcmplx(wr, wi)

	end block

	!print *, "aa = "
	!print "("//to_str(n)//"es16.6)", aa
	!print *, "pq = "
	!print "("//to_str(n)//"es16.6)", pq

	if (.not. present(eigvecs)) return

	! Zero the remaining below-diagonal non-zeros
	call real_schur_to_complex(aa, pq, ca, cq, eps)

	! No eigval swapping happens with LAPACK
	!eigvals = diag(ca)

	!print *, "ca = "
	!print "("//to_str(2*n)//"es19.9)", ca

	! Find the eigenvectors of triangular `ca`
	eigvecs = eye(n)
	do i = 2, n
		eigvecs(1: i-1, i) = invmul(ca(i,i) * eye(i-1) - ca(:i-1, :i-1), ca(:i-1, i), io)
		if (io /= 0) then
			msg = "invmul() failed in eig_lapack()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if
	end do
	!print *, "R eigvecs = "
	!print "("//to_str(n)//"es19.9)", eigvecs

	eigvecs = matmul(cq, eigvecs)

end function eig_lapack

!===============================================================================

function eig_francis_qr(aa, eigvecs, iostat) result(eigvals)
	! Get the eigenvalues of `aa` using the Francis double step QR algorithm.
	! Use eig_lapack() instead because it's much more robust
	use numa__blarg
	use numa__utils
	double precision, intent(inout) :: aa(:,:)
	double complex, allocatable :: eigvals(:)
	double complex, optional, allocatable, intent(out) :: eigvecs(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg

	double complex :: l1, l2
	double complex, allocatable :: ca(:,:), cq(:,:)

	! I still have doubts about robustness here because this tolerance is
	! extremely sensitive
	double precision, parameter :: eps = 1.d-10  ! should eps and iters be args?

	double precision :: rad, s, tr, x, y, z, p2(2,2), p3(3,3), ck, sk, &
		a, b, c, d, det_, v(3), aii, aij, aji, ajj, a11, a12, a21, a22, a32
	double precision, allocatable :: pq(:,:)

	integer, parameter :: iters = 100
	integer :: i, i1, k, n, j, k1, iter, io

	logical, allocatable :: is_real(:)

	if (present(iostat)) iostat = 0

	n = size(aa, 1)

	if (present(eigvecs)) then
		! The matrix `pq` is the product of Q from each step.  Unlike basic QR, we
		! can't initialize pq to eye here.  Instead, we need an initial
		! transformation from the Hessenberg reduction
		call hess(aa, pq)
	else
		call hess(aa)
	end if

	!print *, "aa = "
	!print "("//to_str(n)//"es19.9)", aa

	i = n  ! i indicates the active matrix size
	iter = 0

	! LAPACK has lots of logic to scale to avoid over/under flow which is
	! skipped here and thus not robust

	do while (i > 2)
		iter = iter + 1
		if (iter > iters) then
			msg = "convergence failed in eig_francis_qr()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if
		!print *, "i = ", i

		j = i - 1
		aii = aa(i,i)
		aij = aa(i,j)
		aji = aa(j,i)
		ajj = aa(j,j)

		tr = (aii + ajj) / 1          ! trace
		det_ = ajj * aii - aji * aij  ! determinant

		! Compute first 3 elements of first column of M

		a11 = aa(1,1)
		a21 = aa(2,1)
		a12 = aa(1,2)
		a22 = aa(2,2)
		a32 = aa(3,2)

		x = a11**2 + a12 * a21 - tr * a11 + det_
		y = a21 * (a11 + a22 - tr)
		z = a21 * a32

		x = x
		y = y
		z = z

		do k = 0, i-3

			! Determine the Householder reflector `p3`.  Using a scalar `s`
			! here, and for the Givens rotation below, effectively aids
			! convergence for n >~ 40
			v = [x, y, z]
			s = sum(abs(v))
			v = v / s
			p3 = house(v, io)

			if (io /= 0) then
				!print *, "x, y, z = ", x, y, z
				msg = "house() failed in eig_francis_qr()"
				call PANIC(msg, present(iostat))
				iostat = 2
				return
			end if
			!print *, "p3 = "
			!print "(3es15.5)", p3

			k1 = max(1, k)
			aa(k+1: k+3, k1:) = matmul(p3, aa(k+1: k+3, k1:))

			k1 = min(k+4, i)
			aa(:k1, k+1: k+3) = matmul(aa(:k1, k+1: k+3), p3)

			if (present(eigvecs)) then
				pq(:, k+1: k+3) = matmul(pq(:, k+1: k+3), p3)
			end if

			x = aa(k+2, k+1)
			y = aa(k+3, k+1)
			if (k < i-3) then
				z = aa(k+4, k+1)
			end if

		end do

		! Determine the 2D Givens rotation `p2`

		s = abs(x) + abs(y)
		if (s <= 0) then
			!print *, "x, y = ", x, y
			msg = "matrix is singular in eig_francis_qr()"
			call PANIC(msg, present(iostat))
			iostat = 3
			return
		end if
		x = x / s
		y = y / s

		rad = norm2([x, y])
		ck =  x / rad
		sk = -y / rad
		p2(1,:) = [ck, -sk]
		p2(2,:) = [sk,  ck]

		aa(j:i, i-2:) = matmul(p2, aa(j:i, i-2:))
		aa(:i, i-1:i) = matmul(aa(:i, i-1:i), transpose(p2))

		if (present(eigvecs)) then
			pq(:, i-1:i) = matmul(pq(:, i-1:i), transpose(p2))
		end if

		! Check for convergence

		!print *, "aa(i,j)     = ", aa(i,j)
		!print *, "aa(i-1,j-1) = ", aa(i-1,j-1)

		if (abs(aa(i,j)) < eps * (abs(aa(j,j)) + abs(aa(i,i)))) then
			!print *, "i -= 1"
			aa(i,j) = 0.d0
			i = i - 1
			iter = 0
		else if (abs(aa(i-1, j-1)) < eps * (abs(aa(j-1, j-1)) + abs(aa(j,j)))) then
			!print *, "i -= 2"
			aa(i-1, j-1) = 0.d0
			i = i - 2
			iter = 0
		end if

	end do
	!print *, "aa = "
	!print "("//to_str(n)//"es19.9)", aa
	!print *, "pq = "
	!print "("//to_str(n)//"es19.9)", pq

	! Process 2x2 block along diagonal and get their eigenvalues.  Some
	! may be real, some may be complex
	allocate(is_real(n))
	is_real = .true.
	allocate(eigvals(n))
	do i = 1, n-1
		i1 = i + 1

		! 2x2 block around diagonal
		a = aa(i , i )
		b = aa(i1, i )
		c = aa(i , i1)
		d = aa(i1, i1)

		if (abs(b) <= eps * (abs(a) + abs(d))) cycle

		tr = a + d         ! trace
		det_ = a*d - b*c  ! determinant

		! Eigenvalues of 2x2 block
		!
		! Source:  https://people.math.harvard.edu/~knill/teaching/math21b2004/exhibits/2dmatrices/index.html
		l1 = tr/2 + sqrt(dcmplx(tr**2/4 - det_))
		l2 = tr/2 - sqrt(dcmplx(tr**2/4 - det_))

		!print *, "b  = ", b
		!print *, "l1 = ", l1
		!print *, "l2 = ", l2
		!print *, ""

		is_real(i ) = .false.
		is_real(i1) = .false.

		eigvals(i ) = l1
		eigvals(i1) = l2

	end do

	! Second pass: collect the real eigenvalues not saved in the first pass
	do i = 1, n
		if (.not. is_real(i)) cycle
		eigvals(i) = aa(i,i)
	end do

	!!print *, "eigvals sorted = ", sorted(eigvals)
	!print *, "eigvals = ", eigvals

	!********
	if (.not. present(eigvecs)) return

	! Zero the remaining below-diagonal non-zeros
	call real_schur_to_complex(aa, pq, ca, cq, eps)

	! Some of the eigenvalues get swapped around in real_schur_to_complex()
	eigvals = diag(ca)

	!print *, "ca = "
	!print "("//to_str(2*n)//"es19.9)", ca

	! Find the eigenvectors of triangular `ca`
	eigvecs = eye(n)
	do i = 2, n
		! Looking at LAPACK, it has a special-purpose quasi-triangular solver
		! for this problem.  There are cases for 1x1 blocks, 2x2 blocks, and
		! real or complex eigenvectors, with all complex numbers encoded as
		! pairs of reals, without explicitly doing rsf2csf:
		!
		!     https://netlib.org/lapack/explore-html/d2/d98/group__trevc3_gaee05b7252c5a3b2b935d5a4a6101033d.html
		!
		eigvecs(1: i-1, i) = invmul(ca(i,i) * eye(i-1) - ca(:i-1, :i-1), ca(:i-1, i), io)
		if (io /= 0) then
			msg = "invmul() failed in eig_francis_qr()"
			call PANIC(msg, present(iostat))
			iostat = 4
			return
		end if
	end do
	!print *, "R eigvecs = "
	!print "("//to_str(n)//"es19.9)", eigvecs

	eigvecs = matmul(cq, eigvecs)

	!print *, "eigvecs francis qr = "
	!print "("//to_str(n)//"es19.9)", eigvecs

end function eig_francis_qr

!===============================================================================

subroutine real_schur_to_complex(r, q, cr, cq, eps)
	! Convert real Schur form [r, q] to complex Schur form [cr, cq].  Compare
	! the MATLAB and scipy functions rsf2csf()
	!
	! Matrix `r` is quasi-triangular and `q` is unitary
	use numa__blarg
	use numa__utils
	double precision, intent(in) :: r(:,:), q(:,:)
	double complex, allocatable, intent(out) :: cr(:,:), cq(:,:)
	double precision, intent(in) :: eps
	!********

	double complex :: mu(2), cc, sc, g(2,2)

	double precision :: a, b, c, d, t, det_, rad

	integer :: i

	! Cast to complex
	cr = r
	cq = q

	do i = 2, size(r, 1)
		! Source:
		!
		!     https://github.com/scipy/scipy/blob/3fe8b5088d1b63e7557d26314cf5f40851f46a45/scipy/linalg/_decomp_schur.py#L329

		if (abs(cr(i, i-1)) <= eps * (abs(cr(i-1, i-1)) + abs(cr(i,i)))) then
			cr(i, i-1) = 0
			cycle
		end if

		!print *, "Zeroing at i = ", i
		!print *, "cr(i, i-1) = ", cr(i, i-1)

		! 2x2 block around diagonal
		a = dble(cr(i-1, i-1))
		b = dble(cr(i, i-1))
		c = dble(cr(i-1, i))
		d = dble(cr(i, i))

		t = a + d         ! trace
		det_ = a*d - b*c  ! determinant

		! Eigenvalues of 2x2 block
		mu(1) = t/2 + sqrt(dcmplx(t**2/4 - det_))
		mu(2) = t/2 - sqrt(dcmplx(t**2/4 - det_))
		!print *, "mu = ", mu

		mu(1) = mu(1) - cr(i,i)

		! Maybe this could work as a replacement for dlanv2(), but cc and sc are
		! complex, not real
		rad = norm2c([mu(1), cr(i, i-1)])
		cc = mu(1) / rad
		sc = cr(i, i-1) / rad

		g(1,:) = [conjg(cc), sc]
		g(2,:) = [-sc, cc]

		cr(i-1:i, i-1:) = matmul(g, cr(i-1:i, i-1:))
		cr(:i, i-1:i)   = matmul(cr(:i, i-1:i), transpose(conjg(g)))
		cq(: , i-1:i)   = matmul(cq(: , i-1:i), transpose(conjg(g)))

		cr(i, i-1) = 0

		!print *, "cr = "
		!print "("//to_str(2*n)//"es19.9)", cr
		!print *, ""

	end do

end subroutine real_schur_to_complex

!===============================================================================

function eig_hess_qr_kernel(a, iters, eigvecs, iostat) result(eigvals)
	! Get the real eigenvalues of `a` using `iters` iterations of the Hessenberg
	! QR algorithm
	!
	! To get the eigenvectors, find the kernel (null-space) of A-lambda*I for
	! each eigenvalue lambda

	use numa__utils, only:  sorted
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable :: eigvals(:)
	integer, intent(in) :: iters
	double precision, optional, allocatable, intent(out) :: eigvecs(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	double precision :: h1, h2, rad, givens(2,2), eigval
	double precision, allocatable :: c(:), s(:), a0(:,:), &
		eigvec(:)
	integer :: i, j, k, n, io

	if (present(iostat)) iostat = 0

	a0 = a

	n = size(a, 1)
	allocate(c(n-1), s(n-1))

	call hess(a)

	do i = 1, iters

		do k = 1, n-1
			h1 = a(k, k)
			h2 = a(k+1, k)

			rad = norm2([h1, h2])
			!print *, "rad = ", rad
			c(k) =  h1 / rad
			s(k) = -h2 / rad
			givens(1,:) = [c(k), -s(k)]
			givens(2,:) = [s(k),  c(k)]

			! Ref:  https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
			a(k:k+1, k:) = matmul(givens, a(k:k+1, k:))

		end do
		!print *, "r = "
		!print "(4es15.5)", a

		!! If you stop here, `a` has been overwritten with `r` from its QR
		!! factorization.  You could also compute Q by applying the matmuls in the
		!! loop above to an initial identity matrix, but we don't need Q explicitly
		!! for a Hessenberg QR step
		!stop

		! The rest of the Hessenberg QR step is not part of the QR factorization,
		! rather, it overwrites `a` with R * Q

		! Apply the Givens rotations from the right
		do k = 1, n-1
			givens(1,:) = [ c(k), s(k)]  ! note this is transposed compared to above
			givens(2,:) = [-s(k), c(k)]
			a(1: k+1, k: k+1) = matmul(a(1: k+1, k: k+1), givens)
		end do

	end do
	!print *, "a = "
	!print "(4es15.5)", a

	! Don't sort the eigvals, because that will break the ordering of eigvecs
	!eigvals = sorted(diag(a))

	eigvals = diag(a)

	if (.not. present(eigvecs)) return
	allocate(eigvecs(n,n))

	!print *, "getting eigvecs via kernel ..."

	do i = 1, n
		eigval = eigvals(i)

		! Get eigenvector by finding the null-space of A - eigval * eye
		a = a0
		do j = 1, n
			a(j,j) = a(j,j) - eigval
		end do
		!print *, "a = "
		!print "(4es15.5)", a

		eigvec = lu_kernel(a, io)
		if (io /= 0) then
			msg = "lu_kernel() failed in eig_hess_qr_kernel()"
			call PANIC(msg, present(iostat))
			iostat = 1
			return
		end if

		!********
		!! Find null-space using QR decomposition.  This also works
		!call qr_factor(a, diag_, iostat = io)
		!q = qr_get_q_expl(a, diag_)
		!!print *, "q = "
		!!print "(4es15.5)", q

		!! A - eigval*eye is singular, so the last row of Q will be in its null space
		!! and thus an eigenvector of A
		!!
		!! This probably won't work for eigenvalues with multiplicity.  In
		!! that case, you would need to get the last several rows and set
		!! multiple rows in the output
		!eigvec = q(:,n)
		!!print *, "eigvec = ", eigvec
		!********

		!! Confirm that it's an eigenvec.  All components of this print should be the
		!! same number
		!print *, "eigval = ", matmul(a0, eigvec) / eigvec
		!print *, "eigval = ", matmul(transpose(a0), eigvec) / eigvec

		eigvecs(:,i) = eigvec

	end do

end function eig_hess_qr_kernel

!===============================================================================

end module numa__eig

