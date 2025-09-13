
!> Function callback definitions for testing interpolators, integrators, etc.
module numa__functions

	implicit none

	double precision, parameter, private :: PI = 4 * atan(1.d0)

contains

!===============================================================================

double precision function log_fn(x) result(log_x)

	! The built-in Fortran fn `log` cannot be passed directly as a function
	! callback:
	!
	!      145 |         fx = lagrange_interpolator(xi, log, 11.1d0)
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

!********

double precision function cotd_fn(x)
	double precision, intent(in) :: x
	cotd_fn = 1 / tand(x)
end function cotd_fn

!********

double precision function sin_fn(x)
	double precision, intent(in) :: x
	sin_fn = sin(x)
end function sin_fn

!********

double precision function sqrt_fn(x)
	double precision, intent(in) :: x
	sqrt_fn = sqrt(x)
end function sqrt_fn

!********

double precision function inv_sqrt_fn(x)
	double precision, intent(in) :: x
	inv_sqrt_fn = 1.d0 / sqrt(x)
end function inv_sqrt_fn

!********

double precision function inv_square_fn(x)
	double precision, intent(in) :: x
	inv_square_fn = 1.d0 / x ** 2
end function inv_square_fn

!********

double precision function inv_1px2(x)
	double precision, intent(in) :: x
	inv_1px2 = 1.d0 / (1.d0 + x ** 2)
end function inv_1px2

!********

double precision function inv_1px4(x)
	double precision, intent(in) :: x
	inv_1px4 = 1.d0 / (1.d0 + x ** 4)
end function inv_1px4

!********

double precision function inv_1px1p5(x)
	double precision, intent(in) :: x
	inv_1px1p5 = 1.d0 / (1.d0 + abs(x) ** 1.5d0)
end function inv_1px1p5

!********

double precision function sinx_x(x)
	double precision, intent(in) :: x
	if (x == 0.d0) then
		sinx_x = 1
		return
	end if
	sinx_x = sin(x) / x
end function sinx_x

!********

double precision function exp_nx2(x)
	double precision, intent(in) :: x
	exp_nx2 = exp(-x ** 2)
end function exp_nx2

!********

double precision function fn_example_pg171(x)
	! This is an example integrand from page 171
	double precision, intent(in) :: x
	fn_example_pg171 = &
		5.d0 / (exp(PI) - 2.d0) * exp(2*x) * cos(x)
end function fn_example_pg171

!********

double precision function bessel_3_2p5(tau)
	! This is the integrand for bessel_jn(3, 2.5)
	double precision, intent(in) :: tau
	double precision, parameter :: x = 2.5d0
	integer, parameter :: n = 3

	bessel_3_2p5 = &
		1/PI * cos(n * tau - x * sin(tau))

end function bessel_3_2p5

!********

double precision function rate_fn(x, beta) result(y)
	double precision, intent(in) :: x, beta(:)
	y = beta(1) * x / (beta(2) + x)
end function rate_fn

function rate_dx_fn(x, beta) result(dy)
	! Derivative of rate_fn() wrt x
	double precision, intent(in) :: x, beta(:)
	double precision, allocatable :: dy(:)
	allocate(dy(2))

	dy(1) = x / (beta(2) + x)
	dy(2) = -beta(1) * x / (beta(2) + x)**2

end function rate_dx_fn

double precision function bimodal_fn(x, pk) result(y)
	double precision, intent(in) :: x, pk(:)
	y = pk(1) * exp(-pk(2) * (x - pk(3))**2) &
	  + pk(4) * exp(-pk(5) * (x - pk(6))**2)
end function bimodal_fn

!********

double precision function rosenbrock_banana_beta(dummy, x) result(y)
	! This is a hack for abusing the data-fitting optimizer as a general-purpose
	! optimizer
	double precision, intent(in) :: dummy, x(:)
	double precision, parameter :: a = 1.d0, b = 100.d0

	y = dummy  ! suppress unused arg warning. no other effect

	y = (a - x(1))**2 + b * (x(2) - x(1))**2
	y = -y

end function rosenbrock_banana_beta

double precision function rosenbrock_banana(x) result(y)
	double precision, intent(in) :: x(:)
	double precision, parameter :: a = 1.d0, b = 100.d0
	y = (a - x(1))**2 + b * (x(2) - x(1))**2
end function rosenbrock_banana

double precision function rosenbrock_banana_nd(x) result(y)
	! Note: this is only defined for even size(x)
	double precision, intent(in) :: x(:)
	!********
	integer :: i, n
	n = size(x)
	y = 0
	do i = 1, n/2
		y = y + 100 * (x(2*i-1)**2 - x(2*i))**2 + (x(2*i-1) - 1)**2
	end do
end function rosenbrock_banana_nd

!===============================================================================

double precision function f_nr_ex1(x) result(y)
	! Example for Newton-Raphson:  y = x**2 - a
	double precision, intent(in) :: x
	double precision, parameter :: a = 612
	y = x**2 - a
end function f_nr_ex1
double precision function df_nr_ex1(x) result(dy)
	! Derivative of f_nr_ex1()
	double precision, intent(in) :: x
	dy = 2 * x
end function df_nr_ex1

!********

double precision function f_nr_ex2(x) result(y)
	! Example for Newton-Raphson:  y = cos(x) - x**3
	double precision, intent(in) :: x
	y = cos(x) - x**3
end function f_nr_ex2
double precision function df_nr_ex2(x) result(dy)
	! Derivative of f_nr_ex2()
	double precision, intent(in) :: x
	dy = -sin(x) - 3 * x**2
end function df_nr_ex2

!===============================================================================

end module numa__functions

