
module exercises_m

	use numerical_analysis_m
	use utils_m

	implicit none

contains

!===============================================================================

integer function chapter_2_exercise_2() result(io)

	! 2. Interpolate the function ln x by a quadratic polynomial at x = 10, 11,
	!    12.
	!    (a) Estimate the error commited for x = 11.1 when approximating ln x by
	!        the interpolating polynomial

	double precision :: x, fx
	double precision, allocatable :: xi(:)

	write(*,*) CYAN // "Starting chapter_2_exercise_2()" // COLOR_RESET

	! Support points `xi`
	xi = [10.d0, 11.d0, 12.d0]

	x = 11.1d0

	fx = lagrange_interpolater(xi, log_fn, x)
	print *, "fx   = ", fx
	print *, "f(x) = ", log_fn(x)
	print *, ""

	fx = lagrange_interpolater(xi, exp_fn, x)  ! not from the text book, but let's see what happens
	print *, "fx   = ", fx
	print *, "f(x) = ", exp_fn(x)
	print *, ""

	fx = neville_interpolater(xi, log_fn, x)
	!fx = neville_interpolater_vals(xi, log(xi), x)
	print *, "fx   = ", fx
	print *, "f(x) = ", log_fn(x)
	print *, ""

	! TODO: other exercises for Neville's algorithm, Newton's interpolation
	! formula, etc.

	io = 0

end function chapter_2_exercise_2

!===============================================================================

end module exercises_m

