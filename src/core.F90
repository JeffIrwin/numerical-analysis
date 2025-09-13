
!> Core module with minimal dependencies
module numa__core

	implicit none

	double precision, parameter :: PI = 4 * atan(1.d0)

	double complex, parameter :: IMAG_ = (0.d0, 1.d0)  ! sqrt(-1)

	abstract interface
		! These are function interfaces for passing callbacks

		! Simple real scalar to scalar function, e.g. for interpolation,
		! integration, or scalar root-finding
		function fn_f64_to_f64(x) result(fx)
			double precision, intent(in) :: x
			double precision :: fx
		end function

		! Vector to scalar, e.g. for optimization
		function fn_vec_f64_to_f64(x) result(fx)
			double precision, intent(in) :: x(:)
			double precision :: fx
		end function

		! Vector to vector, e.g. for multi-dimensional root-finding
		function fn_vec_f64_to_vec_f64(x) result(fx)
			double precision, intent(in) :: x(:)
			double precision, allocatable :: fx(:)
		end function

		! Vector to matrix, e.g. for root-finding derivative/Jacobian
		function fn_vec_f64_to_mat_f64(x) result(fx)
			double precision, intent(in) :: x(:)
			double precision, allocatable :: fx(:,:)
		end function

		! Scalar to scalar with additional parameters, e.g. for non-linear least
		! squares in gauss_newton()
		function fn_f64_params_to_f64(x, params) result(fx)
			double precision, intent(in) :: x
			double precision, intent(in) :: params(:)
			double precision :: fx
		end function

		! Scalar to vector with params parameters, e.g. the derivatives of
		! fn_f64_params_to_f64
		function fn_f64_params_to_vec_f64(x, params) result(fx)
			double precision, intent(in) :: x
			double precision, intent(in) :: params(:)
			double precision, allocatable :: fx(:)
		end function

	end interface

end module numa__core

