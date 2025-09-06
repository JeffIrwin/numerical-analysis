
!===============================================================================

#include "panic.F90"

!> Module for high-level linear algebra.  Some other more basic routines can be
!> found in the module numa__blarg
module numa__linalg

	!use numa__core
	use numa__blarg

	implicit none

	interface lu_invmul
		! The built-in `matmul()` works on matrices and vectors similarly
		!
		! This is the subroutine version which modifies its args
		procedure :: lu_invmul_vec
		procedure :: lu_invmul_mat
	end interface lu_invmul

	interface lu_factor
		procedure :: lu_factor_c64
		procedure :: lu_factor_f64
	end interface lu_factor
	interface lu_solve
		procedure :: lu_solve_c64
		procedure :: lu_solve_f64
	end interface lu_solve
	interface lu_kernel
		procedure :: lu_kernel_c64
		procedure :: lu_kernel_f64
	end interface lu_kernel

	interface backsub
		procedure :: backsub_c64
		procedure :: backsub_f64
	end interface backsub

	interface invmul
		! This is the function version which returns the result
		procedure :: invmul_vec_c64
		!procedure :: invmul_mat_c64
		procedure :: invmul_vec_f64
		procedure :: invmul_mat_f64
	end interface invmul

	interface qr_solve
		!procedure :: qr_solve_c64
		procedure :: qr_solve_f64
	end interface qr_solve
	interface qr_factor
		procedure :: qr_factor_c64
		procedure :: qr_factor_f64
	end interface qr_factor
	interface qr_get_q_expl
		procedure :: qr_get_q_expl_c64
		procedure :: qr_get_q_expl_f64
	end interface qr_get_q_expl

	interface qr_mul
		procedure :: qr_mul_mat_c64
		procedure :: qr_mul_mat_f64
		! Needs vec overloads
	end interface qr_mul

	interface qr_mul_transpose
		procedure :: qr_mul_transpose_vec_f64
		procedure :: qr_mul_transpose_mat_f64
		! Needs c64 overloads
	end interface qr_mul_transpose

contains

!===============================================================================

subroutine tridiag_factor(a, iostat)
	! This is the factorization phase of tridiag_invmul(), without pivoting
	!
	! If you are only solving once, the interface for tridiag_invmul() is
	! easier.  However, separate factor and solve phases are better if you need
	! to factor once and then solve multiple times without repeating the
	! factorization
	!
	! Compare the interface of lapack routines dgetrf() (factor) and dgetrs()
	! (solve) vs dgesv() (combined factor-solve).  For example:
	!
	!     https://github.com/JeffIrwin/ribbit/blob/ec868aa3db96258b95909f3101434128fc42428f/src/ribbit.f90#L1722

	use numa__utils
	double precision, intent(inout) :: a(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	integer :: i, n

	if (present(iostat)) iostat = 0

	n = size(a, 2)  ! *not* same as size 1, which is always 3 for this tridiag storage scheme!

	if (size(a, 1) /= 3) then
		msg = "tridiagonal `a` size 1 is not 3 in tridiag_factor()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (a(1,1) /= 0) then
		msg = "`a(1,1)` is not 0 in tridiag_factor()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if
	if (a(3,n) /= 0) then
		msg = "`a(3,n)` is not 0 in tridiag_factor()"
		call PANIC(msg, present(iostat))
		iostat = 3
		return
	end if

	! This implementation is based on wikipedia's "Tridiagonal
	! matrix algorithm" page.  Note that variable names and
	! subscripts differ.

	! l(1) = a(1,1) is 0, u(n) = a(3,n) is 0

	! Could add an `allow_singular` arg like other factor routines
	if (a(2, 2) == 0) then
		msg = "matrix is singular in tridiag_factor()"
		call PANIC(msg, present(iostat))
		iostat = 4
		return
	end if
	a(3,1) = a(3,1) / a(2,1)
	do i = 2, n - 1
		a(2, i) = a(2, i) - a(1, i) * a(3, i-1)
		if (a(2, i) == 0) then
			msg = "matrix is singular in tridiag_factor()"
			call PANIC(msg, present(iostat))
			iostat = 4
			return
		end if
		a(3, i) = a(3, i) / a(2, i)
	end do

end subroutine tridiag_factor

!********

subroutine tridiag_solve(a, bx)
	! This is the solving phase of tridiag_invmul()

	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's

	!********

	integer :: i, n

	n = size(a, 2)  ! *not* same as size 1, which is always 3 for this tridiag storage scheme!

	! Forward sweep
	bx(1)  = bx(1)  / a(2,1)
	do i = 2, n - 1
		bx(i) = (bx(i) - a(1, i) * bx(i-1)) / a(2, i)
	end do
	bx(n) = (bx(n) - a(1,n) * bx(n-1)) / (a(2,n) - a(1,n) * a(3,n-1))

	! Back substitution
	do i = n - 1, 1, -1
		bx(i) = bx(i) - a(3,i) * bx(i+1)
	end do

end subroutine tridiag_solve

!********

subroutine tridiag_invmul(a, bx, iostat)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b
	!
	! Tridiagonal matrix `a` is packed into a 3xn array.  On input, `bx` is the
	! RHS `b`.  On output, `bx` is the solution `x`
	!
	! See exercises.f90 for an example of how to pack the `a` matrix

	use numa__utils
	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's
	integer, optional, intent(out) :: iostat

	! I'm not sure if it's better for performance to use the current
	! representation of `a` or its transpose instead.  Most expressions in
	! factor/solve work on a 3x2 or 2x1 submatrix of `a` currently, so it might
	! be better as-is, but not much of a difference

	!********
	character(len = :), allocatable :: msg
	integer :: io

	if (present(iostat)) iostat = 0
	call tridiag_factor(a, io)
	if (io /= 0) then
		msg = "tridiag_factor() failed in tridiag_invmul()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	call tridiag_solve(a, bx)

end subroutine tridiag_invmul

!===============================================================================
subroutine tridiag_corner_invmul(aa, bx, iostat)

	! Solve an augmented tridiagonal system with off-diagonal opposite corner
	! elements, AKA the cyclic Thomas algorithm
	!
	! The cyclic code on the Tridiagonal_matrix_algorithm wikipedia page is
	! straight up wrong.  Here are better references:
	!
	!   - https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
	!   - https://sachinashanbhag.blogspot.com/2014/03/cyclic-tridiagonal-matrices.html

	use numa__utils
	double precision, intent(inout) :: aa(:,:)
	double precision, intent(inout) :: bx(:)
	integer, optional, intent(out)  :: iostat

	!********
	character(len = :), allocatable :: msg
	double precision :: a1, b1, cn, v1, vn
	double precision, allocatable :: u(:)
	integer :: i, n, io

	if (present(iostat)) iostat = 0

	! Make a modified tridiagonal system with no corners and modified first and
	! last main diagonal elements
	n = size(aa, 2)

	a1 = aa(1, 1)  ! top right elem
	b1 = aa(2, 1)  ! first main diag elem
	cn = aa(3, n)  ! lower left elem

	! Make vectors `u` and `v` for Sherman-Morrison formula
	!
	! We need a scratch vector for u.  Although it is initially sparse, we solve
	! `aa \ u` later which fills it in.  For `v` vector we only need v(1) and
	! v(n)
	u = [(0.d0, i = 1, n)]
	u(1) = -b1;  u(n) = cn
	v1 = 1.d0 ;  vn = -a1 / b1

	aa(1,1) = 0.d0
	aa(3,n) = 0.d0
	aa(2,1) = aa(2,1) - u(1) * v1
	aa(2,n) = aa(2,n) - u(n) * vn

	!print *, "aa = "
	!print "(3es16.6)", aa

	! Solve two pure tridiag problems

	call tridiag_factor(aa, io)
	if (io /= 0) then
		msg = "tridiag_factor() failed in tridiag_corner_invmul()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	call tridiag_solve(aa, bx)
	call tridiag_solve(aa, u)

	!print *, "bx = ", bx
	!print *, "u = ", u

	! Solution to cylcic problem using the Sherman-Morrison formula
	bx = bx - ((v1*bx(1) + vn*bx(n)) / (1.d0 + (v1*u(1) + vn*u(n)))) * u

	!print *, "solution = ", bx

end subroutine tridiag_corner_invmul

!===============================================================================

subroutine banded_factor(a, al, indx, nl, nu, iostat)
	! This is the factorization phase of banded_invmul(), with pivoting
	!
	! Compare the interface of lapack routines dgbtrf() (factor) and dgbtrs()
	! (solve) vs dgbsv() (combined factor-solve)

	use numa__utils

	! Reconsider the order of arguments -- maybe in(out) first and out last.  Or
	! perhaps it's better to leave as-is for consistency with banded_solve()
	! args
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable, intent(out) :: al(:,:)
	integer, allocatable, intent(out) :: indx(:)
	integer, intent(in) :: nl, nu
	integer, optional, intent(out) :: iostat

	!********
	character(len = :), allocatable :: msg
	double precision :: d, dum
	integer :: i, j, k, l, n, mm

	if (present(iostat)) iostat = 0

	!print *, "nl, nu = ", nl, nu

	n = size(a, 2)  ! *not* same as size 1
	mm = nl + 1 + nu

	if (nl + nu + 1 /= size(a, 1)) then
		msg = "bad size of matrix `a` in banded_factor()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if
	if (any(a(1: nl, 1) /= 0)) then
		msg = "`a(1:nl, 1)` is not 0 in tridiag_factor()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if
	!if (any(a(mm-nu+1: mm, n) /= 0)) then
	if (any(a(nl+2: mm, n) /= 0)) then
		msg = "`a(nl+2: mm, n)` is not 0 in tridiag_factor()"
		call PANIC(msg, present(iostat))
		iostat = 3
		return
	end if

	! Left-shift the top nl rows by the number of zeros
	do i = 1, nl
		a(1: mm-nl+i-1, i) = a(1+nl-i+1: mm, i)
		a(mm - nl + i: mm, i) = 0
	end do
	d = 1  ! not used for anything, but it shows even vs odd row pivot permutations

	!print *, "a first loop = "
	!print "(5es18.6)", a

	allocate(indx(n))
	allocate(al(nl, n))
	l = nl
	do k = 1, n
		dum = a(1, k)
		i = k
		if (l < n) l = l + 1

		do j = k+1, l
			if (abs(a(1, j)) > abs(dum)) then
				dum = a(1, j)
				i = j
			end if
		end do

		indx(k) = i
		if (i /= k) then
			d = -d
			a(:, [i, k]) = a(:, [k, i])
		end if

		do i = k+1, l
			dum = a(1,i) / a(1,k)
			al(i-k, k) = dum
			a(1:mm-1, i) = a(2:mm, i) - dum * a(2:mm, k)
			a(mm, i) = 0
		end do
	end do
	!print *, "a = "
	!print "(5es18.6)", a
	!print *, "al = "
	!print "(2es18.6)", al

end subroutine banded_factor

!********

subroutine banded_solve(a, al, indx, nl, nu, bx)

	! This is the solving phase of banded_invmul()

	double precision, intent(in) :: a(:,:)
	double precision, intent(in) :: al(:,:)
	integer, intent(in) :: indx(:)

	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's

	integer, intent(in) :: nl, nu

	!********
	double precision :: dum
	integer :: i, j, k, l, n, mm

	!********
	! Solve stage

	n = size(a, 2)  ! *not* same as size 1
	mm = nl + 1 + nu

	l = nl
	do k = 1, n
		j = indx(k)
		if (j /= k) bx([j, k]) = bx([k, j])
		if (l < n) l = l + 1
		bx(k+1: l) = bx(k+1: l) - al(1: l-k, k) * bx(k)
	end do

	l = 1
	do i = n, 1, -1
		dum = bx(i) - dot_product(a(2:l, i), bx(i+1: l+i-1))
		bx(i) = dum / a(1, i)
		if (l < mm) l = l + 1
	end do

end subroutine banded_solve

!********

!> @brief  Solves a linear system `a*x == b` for `x` with a banded matrix `a`
!>
!> @details
!> Banded matrix `a` is packed into a mostly sparse rectangular array.  For
!> example, this data structure:
!>
!>     nl = 2  ! number of lower bands
!>     nu = 1  ! number of upper bands
!>     n  = 5  ! size of x
!>     allocate(a(nl+nu+1, n))
!>     a(:, 1) = [0, 0, 5, 2]
!>     a(:, 2) = [0, 3, 6, 2]
!>     a(:, 3) = [2, 3, 7, 2]
!>     a(:, 4) = [2, 3, 8, 2]
!>     a(:, 5) = [2, 3, 3, 0]
!>
!> represents this dense matrix in MATLAB/Scilab layout:
!>
!>     a = [
!>         5 2 0 0 0
!>         3 6 2 0 0
!>         2 3 7 2 0
!>         0 2 3 8 2
!>         0 0 2 3 3
!>     ]
!>
!> @param[in,out]  a   Banded matrix on input, factored on output
!> @param[in,out]  bx  RHS `b` on input, `x` on output
!> @param[in]      nl  Number of lower non-zero bands below diagonal
!> @param[in]      nu  Number of upper non-zero bands above diagonal
subroutine banded_invmul(a, bx, nl, nu)
	! Maybe should {a, nl, nu} be a struct?  One of the n*'s could be redundant

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)  ! could be extended to rank-2 for multiple RHS's
	integer, intent(in) :: nl, nu

	!********

	double precision, allocatable :: al(:,:)
	integer, allocatable :: indx(:)

	call banded_factor(a, al, indx, nl, nu)
	call banded_solve (a, al, indx, nl, nu, bx)

end subroutine banded_invmul

!===============================================================================

function invmul_vec_c64(a, b) result(x)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b

	double complex, intent(in) :: a(:,:)
	double complex, intent(in) :: b(:)
	double complex, allocatable :: x(:)

	!********

	double complex, allocatable :: aa(:,:)

	integer :: i
	integer, allocatable :: pivot(:)

	x = b
	aa = a

	! Initialize pivot to identity
	pivot = [(i, i = 1, size(aa,1))]

	! TODO: most calls to lu_factor() should check iostat
	call lu_factor(aa, pivot)

	!print *, "pivot = ", pivot
	!print *, "lu_factor(aa) = "
	!!print "(5es18.6)", aa
	!print "(6es15.5)", transpose(aa)

	call lu_solve(aa, x, pivot)

end function invmul_vec_c64

subroutine lu_invmul_vec(a, bx)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:)

	!********

	integer :: i
	integer, allocatable :: pivot(:)

	! Initialize pivot to identity
	pivot = [(i, i = 1, size(a,1))]

	call lu_factor(a, pivot)

	!print *, "pivot = ", pivot
	!print *, "lu_factor(a) = "
	!!print "(5es18.6)", a
	!print "(6es15.5)", transpose(a)

	call lu_solve(a, bx, pivot)

end subroutine lu_invmul_vec

!********

function invmul_vec_f64(a, b) result(x)
	double precision, intent(in)  :: a(:,:)
	double precision, intent(in)  :: b(:)
	double precision, allocatable :: x(:)
	!********
	double precision, allocatable :: aa(:,:)

	x = b
	aa = a
	call lu_invmul(aa, x)

end function invmul_vec_f64

function invmul_mat_f64(a, b) result(x)
	double precision, intent(in)  :: a(:,:)
	double precision, intent(in)  :: b(:,:)
	double precision, allocatable :: x(:,:)
	!********
	double precision, allocatable :: aa(:,:)

	x = b
	aa = a
	call lu_invmul(aa, x)

end function invmul_mat_f64

!********

subroutine lu_invmul_mat(a, bx)
	! Solve the linear algebra problem for `x`:
	!
	!     a * x = b

	double precision, intent(inout) :: a(:,:)
	double precision, intent(inout) :: bx(:,:)

	!********

	integer :: i, nrhs
	integer, allocatable :: pivot(:)

	! Initialize pivot to identity
	pivot = [(i, i = 1, size(a,1))]

	call lu_factor(a, pivot)

	!print *, "pivot = ", pivot
	!print *, "lu_factor(a) = "
	!!print "(5es18.6)", a
	!print "(6es15.5)", transpose(a)

	nrhs = size(bx, 2)
	do i = 1, nrhs
		call lu_solve(a, bx(:,i), pivot)
	end do

end subroutine lu_invmul_mat

!********

subroutine lu_factor_f64(a, pivot, allow_singular, iostat)

	!! Make pivoting optional (for comparison to other solvers)
	!logical, parameter :: DO_PIVOT = .false.

	use numa__utils
	double precision, intent(inout) :: a(:,:)
	integer, allocatable, intent(inout) :: pivot(:)
	integer, optional, intent(out) :: iostat
	logical, optional, intent(in) :: allow_singular

	!********

	character(len = :), allocatable :: msg
	double precision :: amax
	integer :: i, j, k, n, max_index
	logical :: allow_singular_

	if (present(iostat)) iostat = 0

	allow_singular_ = .false.
	if (present(allow_singular)) allow_singular_ = allow_singular

	n = size(a, 1)

	if (size(a, 2) /= n) then
		msg = "matrix is not square in lu_factor_f64()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	if (.not. allocated(pivot)) then
		! Or we could just allocate and initialize (to identity perm) here
		!
		! Maybe this shouldn't be checked manually at all and we should just let
		! a stacktrace get issued, as with `a`
		msg = "pivot is not allocated in lu_factor_f64()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if
	!print *, "pivot init = ", pivot

	do i = 1, n
		! Find max value in column i
		max_index = i
		amax = abs(a(pivot(i), i))
		!if (DO_PIVOT) then
			do j = i+1, n
				!print *, "aj, ai = ", a(j,i), a(i,i)
				if (abs(a(pivot(j), i)) > amax) then
					max_index = j
					amax = abs(a(pivot(j), i))
				end if
			end do
		!end if

		! Swap rows
		pivot([i, max_index]) = pivot([max_index, i])

		!print *, "swapping ", i, max_index
		!print *, "i = ", i
		!call print_mat(a, "a unpiv = ")
		!call print_mat(a(pivot,:), "a pivot = ")
		!print *, "ap = ", a(pivot(i), i)
		!print *, "ai = ", a(i, i)
		!print *, "pivot = ", pivot
		!print *, ""

		if (.not. allow_singular_ .and. a(pivot(i), i) == 0) then
			msg = "matrix is singular in lu_factor_f64()"
			call PANIC(msg, present(iostat))
			iostat = 3
			return
		end if

		do j = i+1, n
			a    (pivot(j), i) = &
				a(pivot(j), i) / &
				a(pivot(i), i)
			do k = i+1, n
				a    (pivot(j), k) = &
					a(pivot(j), k) - &
					a(pivot(j), i) * &
					a(pivot(i), k)
			end do
		end do

	end do

end subroutine lu_factor_f64

!********

subroutine lu_factor_c64(a, pivot, allow_singular, iostat)

	use numa__utils
	double complex, intent(inout) :: a(:,:)
	integer, allocatable, intent(inout) :: pivot(:)
	logical, optional, intent(in) :: allow_singular
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision :: amax
	integer :: i, j, k, n, max_index
	logical :: allow_singular_

	if (present(iostat)) iostat = 0

	allow_singular_ = .false.
	if (present(allow_singular)) allow_singular_ = allow_singular

	!print *, "pivot = ", pivot

	n = size(a, 1)

	if (size(a, 2) /= n) then
		msg = "matrix is not square in lu_factor_c64()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	if (.not. allocated(pivot)) then
		! Or we could just allocate and initialize (to identity perm) here
		!
		! Maybe this shouldn't be checked manually at all and we should just let
		! a stacktrace get issued, as with `a`
		msg = "pivot is not allocated in lu_factor_c64()"
		call PANIC(msg, present(iostat))
		iostat = 2
		return
	end if

	do i = 1, n
		! Find max value in column i
		max_index = i
		amax = abs(a(pivot(i), i))
		do j = i+1, n
			!if (abs(a(pivot(j),i)) > abs(a(pivot(i),i))) max_index = j
			if (abs(a(pivot(j), i)) > amax) then
				max_index = j
				amax = abs(a(pivot(j), i))
			end if
		end do

		! Swap rows
		pivot([i, max_index]) = pivot([max_index, i])

		if (.not. allow_singular_ .and. a(pivot(i), i) == 0) then
			msg = "matrix is singular in lu_factor_c64()"
			call PANIC(msg, present(iostat))
			iostat = 3
			return
		end if

		do j = i+1, n
			a    (pivot(j), i) = &
				a(pivot(j), i) / &
				a(pivot(i), i)
			do k = i+1, n
				a    (pivot(j), k) = &
					a(pivot(j), k) - &
					a(pivot(j), i) * &
					a(pivot(i), k)
			end do
		end do

	end do

end subroutine lu_factor_c64

!********

subroutine cholesky_factor(a, allow_singular, iostat)
	! Cholesky factorization is only possible for positive definite matrices
	! `a`!  There is no enforcement here

	use numa__utils
	double precision, intent(inout) :: a(:,:)
	logical, optional, intent(in) :: allow_singular
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision :: x, p
	integer :: i, j, k, n
	logical :: allow_singular_

	if (present(iostat)) iostat = 0

	allow_singular_ = .false.
	if (present(allow_singular)) allow_singular_ = allow_singular

	n = size(a, 1)

	if (size(a, 2) /= n) then
		msg = "matrix is not square in cholesky_factor()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	do i = 1, n
	do k = i, n
		x = a(i,k)
		do j = i-1, 1, -1
			x = x - a(k,j) * a(i,j)
		end do
		if (i == k) then
			if (x <= 0.d0 .and. .not. allow_singular_) then
				msg = "matrix is singular in cholesky_factor()"
				call PANIC(msg, present(iostat))
				iostat = 1
				return
			end if
			p = 1.d0 / sqrt(x)
		end if
		a(k,i) = x * p
	end do
	end do

end subroutine cholesky_factor

!********

subroutine qr_factor_gram_schmidt(a, r)
	! Replace `a` with its QR factorization using modified Gram-Schmidt
	!
	! The Householder algorithm in `qr_factor()` is better because it is more
	! numerically stable

	double precision, intent(inout) :: a(:,:)
	double precision, allocatable, intent(out) :: r(:,:)
	!********

	double precision, allocatable :: ai(:)

	integer :: i, j, n

	n = size(a, 1)
	r = zeros(n, n)

	do i = 1, n

		ai = a(:,i)  ! save backup to get r later
		do j = 1, i-1

			a(:,i) = a(:,i) - a(:,j) * &
				dot_product(a(:,j), a(:,i)) / dot_product(a(:,j), a(:,j))

		end do
		a(:,i) = a(:,i) / norm2(a(:,i))

		r(i,i) = dot_product(a(:,i), ai)
		do j = i+1, n
			r(i,j) = dot_product(a(:,i), a(:,j))
		end do

	end do

end subroutine qr_factor_gram_schmidt

!********

function kernel(a, tol, allow_rect, iostat)
	! Get the kernel of `a` using QR factorization
	!
	! Matrix `a` is modified in the process
	!
	! TODO: DRY kernel(), qr_rank(), and qr_factor_f64().  Need to make an
	! optional `pivot` arg, false by default, but true for rank and kernel.  A
	! core routine should provide qr factorization and opt out-arg rank
	use numa__utils
	double precision, intent(inout) :: a(:,:)

	double precision, optional, intent(in) :: tol
	logical, optional, intent(in) :: allow_rect
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: kernel(:,:)

	!********

	character(len = :), allocatable :: msg
	double precision :: s, normx, normj, u1, wa, tol_
	double precision, allocatable :: diag_(:)
	integer :: i, j, k, n, max_index, rank_
	integer, allocatable :: pivot(:)
	logical :: allow_rect_

	allow_rect_ = .false.
	tol_ = 1.d-12
	if (present(iostat)) iostat = 0
	if (present(allow_rect)) allow_rect_ = allow_rect
	if (present(tol)) tol_ = tol

	n = min(size(a,1), size(a,2))

	if (.not. allow_rect_ .and. size(a, 1) /= size(a, 2)) then
		msg = "matrix is not square in kernel()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	allocate(diag_(n))

	! For the right kernel, a transposition is needed.  Otherwise the result is
	! the left kernel
	a = transpose(a)

	! Column pivot, not row pivot unlike LU factor.  Could just transpose
	! everything, might not matter much
	pivot = [(i, i = 1, n)]

	! Ref:  https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	rank_ = 0
	do i = 1, n
		rank_ = rank_ + 1

		max_index = i
		normx = norm2(a(i:, pivot(i)))
		!print *, "normx = ", normx
		do j = i+1, n
			normj = norm2(a(i:, pivot(j)))
			if (normj > normx) then
				max_index = j
				normx = normj
			end if
		end do
		!print *, "normx = ", normx

		! Swap cols
		pivot([i, max_index]) = pivot([max_index, i])

		!if (normx <= 0) then
		if (normx <= tol_) then
			rank_ = rank_ - 1
			cycle
			!exit
		end if

		s = -sign_(a(i, pivot(i)))
		u1 = a(i, pivot(i)) - s * normx
		a(i+1:, pivot(i)) = a(i+1:, pivot(i)) / u1
		a(i, pivot(i)) = s * normx
		diag_(i) = -s * u1 / normx

		do k = i+1, n
			wa = diag_(i) * (dot_product(a(i+1:, pivot(i)), a(i+1:, pivot(k))) + a(i, pivot(k)))
			a(i, pivot(k)) = a(i, pivot(k)) - wa
			a(i+1:, pivot(k)) = a(i+1:, pivot(k)) - a(i+1:, pivot(i)) * wa
		end do

	end do

	!! TODO
	!call qr_factor_f64(a, diag_, allow_singular = .true.)
	!kr = qr_get_q_expl(a, diag_)
	!kr = kr(:, rank_+1:)
	!call print_mat(kr, "kr = ")

	kernel = qr_get_q_expl(a, diag_, pivot)
	kernel = kernel(:, rank_+1:)

end function kernel

!********

integer function qr_rank(a, tol, allow_rect, iostat) result(rank_)
	! Get the rank of `a` using QR factorization with Householder transformations
	!
	! "Rank" here is the meaning in the linear algebra sense, i.e. the dimension
	! of the vector space spanned by the columns of `a`.  This is not to be
	! confused with the Fortran intrinsic sense of "rank", which is a built-in
	! fn being 2 for all matrices
	!
	! TODO: rename just rank_?
	!
	! Matrix `a` is modified in the process
	use numa__utils
	double precision, intent(inout) :: a(:,:)

	double precision, optional, intent(in) :: tol
	logical, optional, intent(in) :: allow_rect
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision :: s, normx, normj, u1, wa, diag_, tol_
	integer :: i, j, k, n, max_index
	integer, allocatable :: pivot(:)
	logical :: allow_rect_

	allow_rect_ = .false.
	tol_ = 1.d-12
	if (present(iostat)) iostat = 0
	if (present(allow_rect)) allow_rect_ = allow_rect
	if (present(tol)) tol_ = tol

	n = min(size(a,1), size(a,2))

	if (.not. allow_rect_ .and. size(a, 1) /= size(a, 2)) then
		msg = "matrix is not square in qr_rank()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	! Column pivot, not row pivot unlike LU factor.  Could just transpose
	! everything, might not matter much
	pivot = [(i, i = 1, n)]

	! Ref:  https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	rank_ = 0
	do i = 1, n
		rank_ = rank_ + 1

		max_index = i
		normx = norm2(a(i:, pivot(i)))
		!print *, "normx = ", normx
		do j = i+1, n
			normj = norm2(a(i:, pivot(j)))
			if (normj > normx) then
				max_index = j
				normx = normj
			end if
		end do
		!print *, "normx = ", normx

		! Swap cols
		pivot([i, max_index]) = pivot([max_index, i])

		!if (normx <= 0) then
		if (normx <= tol_) then
			rank_ = rank_ - 1
			cycle
			!exit
		end if

		s = -sign_(a(i, pivot(i)))
		u1 = a(i, pivot(i)) - s * normx
		a(i+1:, pivot(i)) = a(i+1:, pivot(i)) / u1
		a(i, pivot(i)) = s * normx
		diag_ = -s * u1 / normx

		do k = i+1, n
			wa = diag_ * (dot_product(a(i+1:, pivot(i)), a(i+1:, pivot(k))) + a(i, pivot(k)))
			a(i, pivot(k)) = a(i, pivot(k)) - wa
			a(i+1:, pivot(k)) = a(i+1:, pivot(k)) - a(i+1:, pivot(i)) * wa
		end do

	end do

end function qr_rank

!********

logical function is_full_rank(a, allow_rect, iostat)
	! Determine if matrix `a` is full-rank
	!
	! TODO: tol?
	!
	! Matrix `a` is modified in the process
	use numa__utils
	double precision, intent(inout) :: a(:,:)

	logical, optional, intent(in) :: allow_rect
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision :: s, normx, u1, wa, diag_
	integer :: j, k, n
	logical :: allow_rect_

	allow_rect_ = .false.
	if (present(iostat)) iostat = 0
	if (present(allow_rect)) allow_rect_ = allow_rect

	n = min(size(a,1), size(a,2))

	if (.not. allow_rect_ .and. size(a, 1) /= size(a, 2)) then
		msg = "matrix is not square in is_full_rank()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	! Ref:  https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	is_full_rank = .true.
	do j = 1, n

		normx = norm2(a(j:, j))
		!print *, "normx = ", normx
		if (normx <= 0) then
		!if (normx <= 1.d-10) then
			is_full_rank = .false.
			return
		end if

		s = -sign_(a(j,j))
		u1 = a(j,j) - s * normx
		a(j+1:, j) = a(j+1:, j) / u1
		a(j,j) = s * normx
		diag_ = -s * u1 / normx

		do k = j+1, n
			wa = diag_ * (dot_product(a(j+1:, j), a(j+1:, k)) + a(j,k))
			a(j,k) = a(j,k) - wa
			a(j+1:, k) = a(j+1:, k) - a(j+1:, j) * wa
		end do

	end do

end function is_full_rank

!********

subroutine qr_factor_f64(a, diag_, allow_rect, allow_singular, iostat)
	! Replace `a` with its QR factorization using Householder transformations
	use numa__utils
	double precision, intent(inout) :: a(:,:)

	double precision, allocatable, intent(out) :: diag_(:)
	logical, optional, intent(in) :: allow_rect, allow_singular
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision :: s, normx, u1, wa
	integer :: j, k, n
	logical :: allow_rect_, allow_singular_

	allow_rect_ = .false.
	allow_singular_ = .false.
	if (present(iostat)) iostat = 0
	if (present(allow_rect)) allow_rect_ = allow_rect
	if (present(allow_singular)) allow_singular_ = allow_singular

	n = min(size(a,1), size(a,2))

	allocate(diag_(n))
	diag_ = 0.d0

	if (.not. allow_rect_ .and. size(a, 1) /= size(a, 2)) then
		msg = "matrix is not square in qr_factor_f64()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	! Ref:  https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	do j = 1, n

		normx = norm2(a(j:, j))
		if (normx <= 0 .and. .not. allow_singular_) then
			msg = "matrix is singular in qr_factor_f64()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if

		s = -sign_(a(j,j))
		u1 = a(j,j) - s * normx
		a(j+1:, j) = a(j+1:, j) / u1
		a(j,j) = s * normx
		diag_(j) = -s * u1 / normx

		do k = j+1, n
			wa = diag_(j) * (dot_product(a(j+1:, j), a(j+1:, k)) + a(j,k))
			a(j,k) = a(j,k) - wa
			a(j+1:, k) = a(j+1:, k) - a(j+1:, j) * wa
		end do

	end do
	!print *, "diag_ = "
	!print "(es15.5)", diag_

	! In the end, `R` is formed by zeroing appropriate elements of `a`
	!
	! And then `Q` can be expensively calculated as inv(r) * a (or some similar
	! transpose), or more efficiently using qr_get_q_expl()
	!
	! Results can be checked by verifying that q' * q == eye() or that a = q * r

end subroutine qr_factor_f64

function qr_mul_mat_f64(qr, diag_, x, pivot) result(qx)
	! Implicitly multiply Q * x with a previously computed QR factorization `qr`
	double precision, intent(in) :: qr(:,:), diag_(:), x(:,:)
	double precision, allocatable :: qx(:,:)
	integer, optional, intent(in) :: pivot(:)
	!********
	double precision :: wq

	integer :: j, k, n

	!print *, "qr = "
	!print "(5es15.5)", qr

	n = min(size(qr,1), size(qr,2))
	qx = x  ! could make a subroutine version which replaces x instead

	if (present(pivot)) then
		do j = n, 1, -1
			do k = 1, size(qr, 2)
				wq = diag_(j) * (dot_product(qr(j+1:, pivot(j)), qx(j+1:, k)) + qx(j,k))
				qx(j, k) = qx(j, k) - wq
				qx(j+1:, k) = qx(j+1:, k) - qr(j+1:, pivot(j)) * wq
			end do
		end do

	else
		do j = n, 1, -1
			do k = 1, size(qr, 2)
				wq = diag_(j) * (dot_product(qr(j+1:, j), qx(j+1:, k)) + qx(j,k))
				qx(j, k) = qx(j, k) - wq
				qx(j+1:, k) = qx(j+1:, k) - qr(j+1:, j) * wq
			end do
		end do

	end if

	!print *, "qx = "
	!print "(5es15.5)", qx

end function qr_mul_mat_f64

function qr_mul_transpose_mat_f64(qr, diag_, x) result(qx)
	! Implicitly multiply transpose(Q) * x with a previously computed QR factorization `qr`
	!
	! TODO: add c64 version and c64 test
	double precision, intent(in) :: qr(:,:), diag_(:), x(:,:)
	double precision, allocatable :: qx(:,:)
	!********
	double precision :: wq

	integer :: j, k, n

	!print *, "qr = "
	!print "(5es15.5)", qr

	n = min(size(qr,1), size(qr,2))
	qx = x  ! could make a subroutine version which replaces x instead

	! This loop order is the only difference from qr_mul()
	do j = 1, n
		do k = 1, size(qx, 2)
			wq = diag_(j) * (dot_product(qr(j+1:, j), qx(j+1:, k)) + qx(j,k))
			qx(j, k) = qx(j, k) - wq
			qx(j+1:, k) = qx(j+1:, k) - qr(j+1:, j) * wq
		end do
	end do
	qx = qx(1:n, :)  ! trim

	!print *, "qx = "
	!print "(5es15.5)", qx

end function qr_mul_transpose_mat_f64

function qr_mul_transpose_vec_f64(qr, diag_, x) result(qx)
	! Implicitly multiply transpose(Q) * x with a previously computed QR factorization `qr`
	!
	! TODO: add c64 version
	double precision, intent(in) :: qr(:,:), diag_(:), x(:)
	double precision, allocatable :: qx(:)
	!********
	double precision :: wq

	integer :: j, n

	!print *, "qr = "
	!print "(5es15.5)", qr

	n = min(size(qr,1), size(qr,2))
	qx = x  ! could make a subroutine version which replaces x instead

	! This loop order is the only difference from qr_mul()
	do j = 1, n
		wq = diag_(j) * (dot_product(qr(j+1:, j), qx(j+1:)) + qx(j))
		qx(j) = qx(j) - wq
		qx(j+1:) = qx(j+1:) - qr(j+1:, j) * wq
	end do
	qx = qx(1:n)  ! trim

	!print *, "qx = "
	!print "(5es15.5)", qx

end function qr_mul_transpose_vec_f64

!********

function qr_get_q_expl_f64(qr, diag_, pivot) result(q)
	! Explicitly get the unitary Q matrix from a previously QR-decomposed
	! matrix.  In most cases you will want to avoid using this fn and just
	! implicitly multiply something by Q instead using qr_mul()
	use numa__blarg
	use numa__utils
	double precision, intent(in) :: qr(:,:), diag_(:)
	integer, optional, intent(in) :: pivot(:)
	double precision, allocatable :: q(:,:)

	if (present(pivot)) then
		q = qr_mul(qr, diag_, eye(size(qr,1)), pivot)
	else
		q = qr_mul(qr, diag_, eye(size(qr,1)))
	end if

end function qr_get_q_expl_f64

!********

subroutine qr_factor_c64(a, diag_, iostat)
	! Replace `a` with its QR factorization using Householder transformations
	use numa__blarg
	use numa__utils
	double complex, intent(inout) :: a(:,:)

	double complex, allocatable, intent(out) :: diag_(:)
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double complex :: wa
	integer :: j, k, n, io

	if (present(iostat)) iostat = 0

	n = size(a, 1)

	if (size(a, 2) /= n) then
		! This could have an allow_rect override like qr_factor_f64()
		msg = "matrix is not square in qr_factor_c64()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	allocate(diag_(n))
	diag_ = 0.d0

	! Ref:  https://www.cs.cornell.edu/~bindel/class/cs6210-f09/lec18.pdf
	do j = 1, n

		call house_c64(a(j,j), a(j+1:, j), diag_(j), io)
		if (io /= 0) then
			msg = "matrix is singular in qr_factor_c64()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if

		do k = j+1, n
			! Note diag_ is conjugated here, unlike in qr_mul_c64()
			wa = conjg(diag_(j)) * (dot_product(a(j+1:, j), a(j+1:, k)) + a(j,k))
			a(j,k) = a(j,k) - wa
			a(j+1:, k) = a(j+1:, k) - a(j+1:, j) * wa
		end do

	end do
	!print *, "diag_ = "
	!print "(es15.5)", diag_

end subroutine qr_factor_c64

function qr_mul_mat_c64(qr, diag_, x) result(qx)
	! Implicitly multiply Q * x with a previously computed QR factorization `qr`
	double complex, intent(in) :: qr(:,:), diag_(:), x(:,:)
	double complex, allocatable :: qx(:,:)
	!********
	double complex :: wq
	!double complex, allocatable :: w(:)
	integer :: j, k, n

	!print *, "qr = "
	!print "(5es15.5)", qr

	n = size(qr, 1)
	qx = x

	! To multiply by transpose(Q) instead, just loop from 1 up to n instead
	do j = n, 1, -1
		! Note that diag_ is *not* conjugated here, unlike in qr_factor_c64()
		do k = 1, n
			wq = diag_(j) * (dot_product(qr(j+1:, j), qx(j+1:, k)) + qx(j,k))
			qx(j, k) = qx(j, k) - wq
			qx(j+1:, k) = qx(j+1:, k) - qr(j+1:, j) * wq
		end do
	end do
	!print *, "qx = "
	!print "(5es15.5)", qx

end function qr_mul_mat_c64

!********

function qr_get_q_expl_c64(qr, diag_) result(q)
	! Explicitly get the unitary Q matrix from a previously QR-decomposed
	! matrix.  In most cases you will want to avoid using this fn and just
	! implicitly multiply something by Q instead using qr_mul()
	use numa__blarg
	use numa__utils
	double complex, intent(in) :: qr(:,:), diag_(:)
	double complex, allocatable :: q(:,:)

	q = qr_mul(qr, diag_, dcmplx(eye(size(qr,1))))

end function qr_get_q_expl_c64

!********

subroutine cholesky_solve(a, bx)

	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)

	integer :: i, j, n

	n = size(a, 1)
	!print *, "bx init = ", bx

	! This is the same as lu_solve(), except:
	! - The pivot is the identity
	! - One of the `a` terms is transposed, because U == transpose(L)
	! - Forward substitution has an extra division step because U has non-zero diagonals

	! Forward substitution
	do i = 1, n
		do j = 1, i-1
			bx(i) = bx(i) - a(i, j) * bx(j)
		end do
		bx(i) = bx(i) / a(i, i)
	end do
	!print *, "bx forward = ", bx

	! Back substitution
	do i = n, 1, -1
		do j = i+1, n
			! Note `a` indices are transposed because its symmetric
			bx(i) = bx(i) - a(j, i) * bx(j)
		end do
		bx(i) = bx(i) / a (i, i)
	end do
	!print *, "bx back = ", bx

	! No pivoting for Cholesky

end subroutine cholesky_solve

!********

subroutine lu_solve_c64(a, bx, pivot)

	double complex, intent(in) :: a(:,:)
	double complex, intent(inout) :: bx(:)
	integer, intent(in) :: pivot(:)

	integer :: i, j, n

	!print *, "pivot lu_solve_c64 = ", pivot

	n = size(a, 1)
	!print *, "bx init = ", bx

	! Forward substitution
	do i = 1, n
		do j = 1, i-1
			bx    (pivot(i)) = &
				bx(pivot(i)) - &
				a (pivot(i), j) * &
				bx(pivot(j))
		end do
	end do
	!print *, "bx forward = ", bx

	call backsub(a, bx, pivot)

	! Unpivot
	bx = bx(pivot)
	!print *, "bx unpivot = ", bx

end subroutine lu_solve_c64

!********

subroutine lu_solve_f64(a, bx, pivot)

	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)
	integer, intent(in) :: pivot(:)

	integer :: i, j, n

	!print *, "pivot lu_solve_f64 = ", pivot

	n = size(a, 1)
	!print *, "bx init = ", bx

	! Forward substitution
	do i = 1, n
		do j = 1, i-1
			bx    (pivot(i)) = &
				bx(pivot(i)) - &
				a (pivot(i), j) * &
				bx(pivot(j))
		end do
	end do
	!print *, "bx forward = ", bx

	call backsub(a, bx, pivot)

	! Unpivot
	bx = bx(pivot)
	!print *, "bx unpivot = ", bx

end subroutine lu_solve_f64

!===============================================================================

function inv(a)
	! Return the inverse of matrix `a` without modifying `a`
	double precision, allocatable, intent(in) :: a(:,:)
	double precision, allocatable :: inv(:,:)

	!double precision, allocatable :: copy_(:,:)

	!! Use LU decomposition
	!copy_ = a
	!inv = eye(size(a,1))
	!call invmul(copy_, inv)

	! Use Gauss-Jordan.  I suspect it does less copying than my LU inverse
	! implementation.  It's around 40% faster, although both methods are O(n**3)
	inv = a
	call gauss_jordan(inv)

end function inv

!===============================================================================

subroutine invert(a)
	! Replace `a` with its inverse
	double precision, allocatable, intent(inout) :: a(:,:)
	!double precision, allocatable :: copy_(:,:)

	!********
	! Use Gauss-Jordan
	call gauss_jordan(a)

	!********
	! Use LU

	!! One move and one copy
	!call move_alloc(a, copy_)
	!a = eye(size(copy_, 1))
	!call invmul(copy_, a)

	!!! Two copies
	!!copy_ = a
	!!a = eye(size(a,1))
	!!call invmul(copy_, a)

	!!! Also two copies
	!!copy_ = eye(size(a,1))
	!!call invmul(a, copy_)
	!!a = copy_

	!!! Break tests
	!!a(1,1) = 1.01d0 * a(1,1)

end subroutine invert

!===============================================================================

subroutine gauss_jordan(a, iostat)
	! Use the Gauss-Jordan method to replace matrix `a` with its inverse
	use numa__utils
	double precision, allocatable, intent(inout) :: a(:,:)
	integer, optional, intent(out) :: iostat

	!********

	character(len = :), allocatable :: msg
	double precision :: max_, hr
	integer, allocatable :: p(:)
	integer :: i, j, k, n, r, r1(1)

	if (present(iostat)) iostat = 0

	n = size(a, 1)
	!print *, "n = ", n
	p = [(i, i = 1, n)]  ! pivot

	if (size(a, 2) /= n) then
		msg = "matrix is not square in gauss_jordan()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	do j = 1, n

		if (j < n) then
			r1 = maxloc(abs(a(j+1: n, j))) + j
		else
			r1 = n
		end if

		!print *, "a : = ", a(j+1: n, j)
		!print *, "maxloc = ", maxloc(abs(a(j+1: n, j)))
		!print *, "r1 = ", r1
		!print *, ""

		r = r1(1)
		max_ = abs(a(r, j))
		if (max_ <= 0) then
			msg = "matrix is singular in gauss_jordan()"
			call PANIC(msg, present(iostat))
			iostat = 2
			return
		end if

		if (r > j) then
			! Swap rows
			a([j, r], :) = a([r, j], :)
			p([j, r]) = p([r, j])
		end if

		! Transform
		hr = 1.d0 / a(j,j)
		a(:,j) = a(:,j) * hr
		a(j,j) = hr
		do k = 1, n
			if (k == j) cycle
			do i = 1, n
				if (i == j) cycle

				a(i,k) = a(i,k) - a(i,j) * a(j,k)
			end do
			a(j,k) = -hr * a(j,k)
		end do

	end do

	! Swap columns
	a(:,p) = a

	!print *, "In Gauss-Jordan:"
	!print *, "a = "
	!print "(5es18.6)", a

end subroutine gauss_jordan

!===============================================================================

function matmul_triu_ge(upper, ge, iostat) result(res)
	! Multiply an upper-triangular matrix by a general dense matrix
	use numa__utils
	double precision, intent(in) :: upper(:,:), ge(:,:)
	double precision, allocatable :: res(:,:)
	integer, optional, intent(out) :: iostat
	!********

	character(len = :), allocatable :: msg
	integer :: i, j, k, ni, nj, nk

	if (present(iostat)) iostat = 0

	ni = size(upper, 1)
	nj = size(upper, 2)
	nk = size(ge, 2)

	if (nj /= size(ge, 1)) then
		msg = "inner dimensions do not agree in matmul_triu_ge()"
		call PANIC(msg, present(iostat))
		iostat = 1
		return
	end if

	res = zeros(ni, nk)
	do k = 1, nk
	do j = 1, nj
	do i = 1, j  ! this loop is half the extent of general matmul

		! Is it possible to overwrite the result to upper?  I don't see how,
		! because it appears on both the LHS and RHS of this expression
		res(i, k) = res(i, k) + upper(i, j) * ge(j, k)

	end do
	end do
	end do

end function matmul_triu_ge

!===============================================================================

function lu_kernel_f64(a) result(kernel)
	! Find the kernel or null-space of `a` using LU decomposition.  This only
	! returns 1 vector of the kernel, not the whole kernel basis
	use numa__utils
	double precision, intent(inout) :: a(:,:)
	double precision, allocatable :: kernel(:)

	integer :: i, n
	integer, allocatable :: pivot(:)

	n = size(a, 1)
	kernel = zeros(n)

	! The kernel of `a` is the same as `u`

	pivot = [(i, i = 1, size(a,1))]
	call lu_factor(a, pivot, allow_singular = .true.)

	! Partial backsub
	kernel(n) = 1.d0
	do i = n-1, 1, -1
		kernel(i) = &
			-dot_product(kernel(i+1:), a(pivot(i), i+1:n)) / a(pivot(i),i)
	end do

end function lu_kernel_f64

!===============================================================================

function lu_kernel_c64(a) result(kernel)
	! Find the kernel or null-space of `a` using LU decomposition
	use numa__utils
	double complex, intent(inout) :: a(:,:)
	double complex, allocatable :: kernel(:)

	integer :: i, n
	integer, allocatable :: pivot(:)

	n = size(a, 1)
	kernel = zeros(n)

	! The kernel of `a` is the same as `u`

	pivot = [(i, i = 1, size(a,1))]
	call lu_factor(a, pivot, allow_singular = .true.)

	! Partial backsub
	kernel(n) = 1.d0
	do i = n-1, 1, -1
		kernel(i) = &
			-dot_product(conjg(kernel(i+1:)), a(pivot(i), i+1:n)) / a(pivot(i), i)
	end do

end function lu_kernel_c64

!===============================================================================

function qr_solve_f64(a, b, allow_rect, iostat) result(x)
	! This can be a regular linear system solver, or with allow_rect, a
	! least-squares solver
	use numa__utils
	double precision, intent(inout) :: a(:,:)
	double precision, intent(in) :: b(:)
	logical, optional, intent(in) :: allow_rect
	integer, optional, intent(out) :: iostat
	double precision, allocatable :: x(:)
	!********

	character(len = :), allocatable :: msg
	double precision, allocatable :: diag_(:)
	integer :: io
	logical :: allow_rect_

	if (present(iostat)) iostat = 0

	allow_rect_ = .false.
	if (present(allow_rect)) allow_rect_ = allow_rect

	call qr_factor(a, diag_, allow_rect_, iostat = io)
	if (io /= 0) then
		msg = "qr_factor() failed in qr_solve_f64()"
		call PANIC(msg, present(iostat))
		iostat = 4
		return
	end if
	!print *, "qr(a) = "
	!print "("//to_str(n+1)//"es16.6)", transpose(a)

	x = qr_mul_transpose(a, diag_, b)
	!print *, "qty = ", x
	call backsub(a, x)

end function qr_solve_f64

!===============================================================================

subroutine backsub_f64(a, bx, pivot)
	! Perform the back-substitution solve phase without pivoting.  This works as
	! a full solve if matrix `a` is upper (right) triangular
	double precision, intent(in) :: a(:,:)
	double precision, intent(inout) :: bx(:)
	integer, optional, intent(in) :: pivot(:)
	!********
	integer :: i, j, n

	n = min(size(a,1), size(a,2))
	if (present(pivot)) then

		do i = n, 1, -1
			do j = i+1, n
				bx    (pivot(i)) = &
					bx(pivot(i)) - &
					a (pivot(i), j) * &
					bx(pivot(j))
			end do
			bx    (pivot(i)) = &
				bx(pivot(i)) / &
				a (pivot(i), i)
		end do

	else

		do i = n, 1, -1
			do j = i+1, n
				bx(i) = bx(i) - a(i, j) * bx(j)
			end do
			bx(i) = bx(i) / a(i, i)
		end do

	end if
	!print *, "bx back = ", bx

end subroutine backsub_f64

!===============================================================================

subroutine backsub_c64(a, bx, pivot)
	! Perform the back-substitution solve phase without pivoting.  This works as
	! a full solve if matrix `a` is upper (right) triangular
	double complex, intent(in) :: a(:,:)
	double complex, intent(inout) :: bx(:)
	integer, optional, intent(in) :: pivot(:)
	!********
	integer :: i, j, n

	n = min(size(a,1), size(a,2))
	if (present(pivot)) then

		do i = n, 1, -1
			do j = i+1, n
				bx    (pivot(i)) = &
					bx(pivot(i)) - &
					a (pivot(i), j) * &
					bx(pivot(j))
			end do
			bx    (pivot(i)) = &
				bx(pivot(i)) / &
				a (pivot(i), i)
		end do

	else

		do i = n, 1, -1
			do j = i+1, n
				bx(i) = bx(i) - a(i, j) * bx(j)
			end do
			bx(i) = bx(i) / a(i, i)
		end do

	end if
	!print *, "bx back = ", bx

end subroutine backsub_c64

!===============================================================================

end module numa__linalg

