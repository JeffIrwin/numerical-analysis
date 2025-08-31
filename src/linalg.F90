
!===============================================================================

#include "panic.F90"

!> Main public-facing module for numerical analysis
module numa__linalg

	!use numa__core

	implicit none

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

end module numa__linalg

