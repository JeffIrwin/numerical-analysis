
module exercises_m

	use numerical_analysis_m
	use utils_m

	implicit none

contains

!===============================================================================

integer function assert(condition)
	! Return 0 if true or 1 if false
	logical, intent(in) :: condition
	if (condition) then
		assert = 0
	else
		assert = 1
	end if
end function assert

!===============================================================================

subroutine test(val, expect, tol, nfail, msg)
	double precision, intent(in) :: val, expect, tol
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg

	!********
	integer :: io

	io = assert(abs(val - expect) < tol)
	nfail = nfail + io
	if (io /= 0) write(*,*) ERROR // "assert failed for value " // &
		to_str(val) // " in "// msg

end subroutine test

!********

subroutine tests(vals, expects, tol, nfail, msg)
	double precision, intent(in) :: vals(:), expects(:), tol
	integer, intent(inout) :: nfail
	character(len = *), intent(in) :: msg
	!********
	integer :: i
	do i = 1, size(vals, 1)
		call test(vals(i), expects(i), tol, nfail, msg)
	end do
end subroutine tests

!===============================================================================

integer function chapter_2_exercise_2() result(nfail)

	! 2. Interpolate the function ln x by a quadratic polynomial at x = 10, 11,
	!    12.
	!    (a) Estimate the error commited for x = 11.1 when approximating ln x by
	!        the interpolating polynomial

	double precision :: x, fx
	double precision, parameter :: tol = 1.d-3
	!double precision, parameter :: tol = 4.d-5
	double precision, allocatable :: xi(:), xs(:), fxs(:)

	write(*,*) CYAN // "Starting chapter_2_exercise_2()" // COLOR_RESET

	nfail = 0

	! Support points `xi`
	xi = [10.d0, 11.d0, 12.d0]

	x = 11.1d0

	! The acceptable tolerance depends on xi, x, and the fn being interpolated

	fx = lagrange_interpolater(xi, log_fn, x)
	!nfail = nfail + assert(abs(fx - log_fn(x)) < tol)
	call test(fx, log_fn(x), tol, nfail, "lagrange_interpolater()")
	print *, "fx   = ", fx, " (Lagrange)"
	print *, "f(x) = ", log_fn(x)
	print *, ""

	!fx = lagrange_interpolater(xi, exp_fn, x)  ! not from the text book, but let's see what happens
	!print *, "fx   = ", fx
	!print *, "f(x) = ", exp_fn(x)
	!print *, ""

	fx = neville_interpolater(xi, log_fn, x)
	!fx = neville_interpolater_vals(xi, log(xi), x)
	call test(fx, log_fn(x), tol, nfail, "neville_interpolater()")
	print *, "fx   = ", fx, " (Neville)"
	print *, "f(x) = ", log_fn(x)
	print *, ""

	xs = [11.1d0, 11.2d0]
	!xs = [11.1d0, 11.2d0, 11.3d0, 11.4d0]
	fxs = newton_interpolater(xi, log_fn, xs)
	!fxs = newton_interpolater_vals(xi, log(xi), xs)
	call tests(fxs, log(xs), tol, nfail, "newton_interpolater()")
	print *, "fxs   = ", fxs, " (Newton)"
	print *, "f(xs) = ", log(xs)
	print *, ""

end function chapter_2_exercise_2

!===============================================================================

integer function chapter_2_example_2p2p4() result(nfail)

	! Compare rational and polynomial interpolation, from the example in section
	! 2.2.4 on page 73
	!
	! For the function f(x) = cot(x) th evalues at [1: 5] degrees have been
	! tabulated.  Determine an approximate value for cot(2.5 degrees) using both
	! polynomial and rational interpolation

	double precision :: x, fx
	double precision, parameter :: tol = 1.d-3
	double precision, allocatable :: xi(:)

	write(*,*) CYAN // "Starting chapter_2_example_2p2p4()" // COLOR_RESET

	nfail = 0

	! Support points `xi`
	xi = [1, 2, 3, 4, 5]

	x = 2.5d0

	! The acceptable tolerance depends on xi, x, and the fn being interpolated

	! TODO: cotand (and cotan) is a GNU extension.  Write a cotd_fn() wrapper
	! based on intrinsic tand()

	! The expected value for polynomial interpolation, 22.6352, is quite far
	! from the actual cot value because cot has a singularity at x == 0.
	! Rational interpolation works better near singularities

	!fx = neville_interpolater(xi, log_fn, x)
	fx = neville_interpolater_vals(xi, cotand(xi), x)
	!call test(fx, cotand(x), tol, nfail, "neville_interpolater()")
	call test(fx, 22.6352d0, tol, nfail, "neville_interpolater()")
	print *, "fx   = ", fx, " (Neville polynomial)"
	print *, "f(x) = ", cotand(x)
	print *, ""

	!fx = neville_rational_interpolater(xi, log_fn, x)
	fx = neville_rational_interpolater_vals(xi, cotand(xi), x)
	call test(fx, cotand(x), tol, nfail, "neville_rational_interpolater()")
	!call test(fx, 22.6352d0, tol, nfail, "neville_rational_interpolater()")
	print *, "fx   = ", fx, " (Neville rational)"
	print *, "f(x) = ", cotand(x)
	print *, ""

end function chapter_2_example_2p2p4

!===============================================================================

end module exercises_m

