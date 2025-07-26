
module numa__functions

	implicit none

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

double precision function sinx_x(x)
	double precision, intent(in) :: x
	if (x == 0.d0) then
		sinx_x = 1
		return
	end if
	sinx_x = sin(x) / x
end function sinx_x

!********

double precision function fn_example_pg171(x)
	! This is an example integrand from page 171
	double precision, intent(in) :: x
	double precision, parameter :: PI = 4 * atan(1.d0)
	fn_example_pg171 = &
		5.d0 / (exp(PI) - 2.d0) * exp(2*x) * cos(x)
end function fn_example_pg171

!===============================================================================

end module numa__functions

