
module numerical_analysis_m

	implicit none

	double precision, parameter :: PI = 4 * atan(1.d0)

	abstract interface
		function fn_f64_to_f64(x) result(fx)
			double precision, intent(in) :: x
			double precision :: fx
		end function
	end interface

contains

integer function chapter_2_exercise_2_wet() result(io)

	! 2. Interpolate the function ln x by a quadratic polynomial at x = 10, 11,
	!    12.
	!    (a) Estimate the error commited for x = 11.1 when approximating ln x by
	!        the interpolating polynomial

	double precision :: x, fx, prod
	double precision, allocatable :: xi(:), fi(:)

	integer :: i, k, n, n1

	! Support points `xi`
	xi = [10.d0, 11.d0, 12.d0]

	! Degree of polynomial + 1 `n1` and the actual degree `n`
	n1 = size(xi, 1)
	n = n1 - 1

	! Function values `fi` at those support points
	!
	! TODO: refactor this to take a general fn pointer/callback
	allocate(fi(n1))
	do i = 1, n1
		fi(i) = log(xi(i))
	end do

	print *, "xi = ", xi
	print *, "fi = ", fi

	! Interpolate at f(x) at `x`
	x = 11.1d0

	! Lagrange interpolation
	fx = 0.d0
	do i = 1, n1
		prod = 1.d0
		do k = 1, n1
			if (i == k) cycle
			prod = prod * (x - xi(k)) / (xi(i) - xi(k))
		end do
		fx = fx + fi(i) * prod
	end do

	print *, "fx = ", fx

	io = 0

end function chapter_2_exercise_2_wet

!===============================================================================
!double precision function lagrange_interpolater(xi, fi, f, x) result(fx)
double precision function lagrange_interpolater(xi, f, x) result(fx)

	! Given support points `xi` and function `f`, interpolate f at `x` using
	! Lagrange interpolation
	!
	! "Lagrange's formula is, in general, not as suitable for actual
	! calculations as some other methods to be described below ... ."

	double precision, intent(in) :: xi(:)
	procedure(fn_f64_to_f64) :: f
	double precision, intent(in) :: x

	!********

	double precision :: prod
	double precision, allocatable :: fi(:)

	integer :: i, k, n1

	! Degree of polynomial + 1
	n1 = size(xi, 1)

	! Function values `fi` at those support points
	allocate(fi(n1))
	do i = 1, n1
		!fi(i) = log(xi(i))
		fi(i) = f(xi(i))
	end do

	print *, "xi = ", xi
	print *, "fi = ", fi

	! Interpolate at f(x) at `x`

	! Lagrange interpolation
	fx = 0.d0
	do i = 1, n1
		prod = 1.d0
		do k = 1, n1
			if (i == k) cycle
			prod = prod * (x - xi(k)) / (xi(i) - xi(k))
		end do
		fx = fx + fi(i) * prod
	end do

	print *, "fx     = ", fx
	print *, "actual = ", f(x)
	print *, ""

end function lagrange_interpolater

!===============================================================================

double precision function log_fn(x) result(log_x)

	! The built-in Fortran fn `log` cannot be passed directly as a function
	! callback:
	!
	!      145 |         fx = lagrange_interpolater(xi, log, 11.1d0)
	!          |                                          1
	!    Error: Symbol ‘log’ at (1) has no IMPLICIT type
	!
	! So we define this wrapper

	double precision, intent(in) :: x

	log_x = log(x)

end function log_fn

!********

double precision function exp_fn(x) result(exp_x)
	double precision, intent(in) :: x
	exp_x = exp(x)
end function exp_fn

!===============================================================================

integer function chapter_2_exercise_2() result(io)

	! 2. Interpolate the function ln x by a quadratic polynomial at x = 10, 11,
	!    12.
	!    (a) Estimate the error commited for x = 11.1 when approximating ln x by
	!        the interpolating polynomial

	double precision :: x, fx, prod
	double precision, allocatable :: xi(:), fi(:)

	! Support points `xi`
	xi = [10.d0, 11.d0, 12.d0]

	!fx = lagrange_interpolater(xi, log, x)
	fx = lagrange_interpolater(xi, log_fn, 11.1d0)

	fx = lagrange_interpolater(xi, exp_fn, 11.1d0)

	io = 0

end function chapter_2_exercise_2

!===============================================================================

end module numerical_analysis_m

